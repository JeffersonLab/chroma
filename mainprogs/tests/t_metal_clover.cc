// Real-gauge regression for the two native SystemSolver factories.
#include "actions/ferm/fermacts/eoprec_clover_fermact_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/invert/metal_solvers/syssolver_linop_metal_clover.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "util/gauge/reunit.h"
#ifdef BUILD_SB
#include "util/ferm/mgproton.h"
#endif
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

using namespace Chroma;
namespace {
  using Q = multi1d<LatticeColorMatrix>;
  using T = LatticeFermion;
  uint64_t hash(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
  }
  void source(T &x, int tag) {
    auto dims = Layout::lattSize();
    size_t volume = 1;
    for (int d = 0; d < 4; ++d)
      volume *= dims[d];
    multi1d<int> coord(4);
    double scale = 1 / std::sqrt(double(volume) * 12);
    for (size_t i = 0; i < volume; ++i) {
      size_t q = i;
      for (int mu = 0; mu < 4; ++mu) {
        coord[mu] = q % dims[mu];
        q /= dims[mu];
      }
      int site = Layout::linearSiteIndex(coord);
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto v = hash(i * 12 + c * 4 + s + uint64_t(tag) * 771337);
          auto &z = x.elem(site).elem(s).elem(c);
          z.real() = (v & 1) ? scale : -scale;
          z.imag() = 0;
        }
    }
  }
  void inverter(XMLBufferWriter &xml, const CloverFermActParams &p,
                const multi1d<int> &blocks) {
    push(xml, "InvertParam");
    write(xml, "RsdTarget", Real(1e-12));
    write(xml, "MaxIter", 65536);
    write(xml, "Blocks", blocks);
    write(xml, "SecondFilter", 16);
    write(xml, "BottomIterations", 4);
    push(xml, "CloverParams");
    write(xml, "Mass", p.Mass);
    write(xml, "u0", p.u0);
    write(xml, "AnisoParam", p.anisoParam);
    if (p.anisoParam.anisoP) {
      write(xml, "clovCoeffR", p.clovCoeffR);
      write(xml, "clovCoeffT", p.clovCoeffT);
    } else
      write(xml, "clovCoeff", p.clovCoeffR);
    pop(xml);
    pop(xml);
  }
  double relative(const LinearOperator<T> &A, const T &x, const T &b) {
    T r = zero;
    A(r, x, PLUS);
    r[A.subset()] = b - r;
    double denominator = toDouble(norm2(b, A.subset()));
    return std::sqrt(toDouble(norm2(r, A.subset())) /
                     (denominator > 0 ? denominator : 1));
  }
} // namespace
int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: t_metal_clover GAUGE_LIME ACTION_XML "
                 "BLOCKS_COMMA_SEPARATED\n";
    return 2;
  }
  std::string gauge_path = argv[1], action_path = argv[2], block_arg = argv[3];
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    static_assert(BASE_PRECISION == 64, "requires FP64 Chroma");
    XMLReader metadata(action_path);
    multi1d<int> dims;
    read(metadata, "/DBMetaData/lattSize", dims);
    Layout::setLattSize(dims);
    Layout::create();
    if (Layout::numNodes() != 1)
      throw std::runtime_error("run the Metal harness on one rank");
    const std::string path = "/DBMetaData/Params/Propagator/FermionAction";
    CloverFermActParams params(metadata, path);
    Q u(4);
    XMLReader file_xml, record_xml;
    QDPFileReader reader(file_xml, gauge_path, QDPIO_SERIAL);
    read(reader, record_xml, u);
    if (reader.bad())
      throw std::runtime_error("LIME gauge/checksum failure");
    reader.close();
    for (int mu = 0; mu < 4; ++mu)
      reunit(u[mu]);
    CreateFermStateEnv::registerAll();
    LinOpSysSolverEnv::registerAll();
    auto create = CreateFermStateEnv::reader(metadata, path);
    Handle<FermState<T, Q, Q>> state((*create)(u));

    Handle<LinearOperator<T>> full(new UnprecCloverLinOp(state, params));
    auto *eo = new EvenOddPrecCloverLinOp(state, params);
    Handle<LinearOperator<T>> odd(eo);
    multi1d<int> blocks(4);
    std::stringstream input(block_arg);
    for (int mu = 0; mu < 4; ++mu) {
      char comma;
      if (mu && (!(input >> comma) || comma != ','))
        throw std::runtime_error("block syntax");
      if (!(input >> blocks[mu]) || blocks[mu] < 1)
        throw std::runtime_error("block extent");
    }
    XMLBufferWriter inv;
    inverter(inv, params, blocks);
    XMLReader xml(inv);
    std::cout << std::setprecision(17);
    T b = zero, reference = zero;
    source(b, 20261112);
    for (bool unpreconditioned : {true, false}) {
      Handle<LinearOperator<T>> A = unpreconditioned ? full : odd;
      T rhs = b;
      if (!unpreconditioned) {
        T inv_even = zero, hop = zero;
        eo->evenEvenInvLinOp(inv_even, b, PLUS);
        eo->oddEvenLinOp(hop, inv_even, PLUS);
        rhs[rb[1]] = b - hop;
      }
      for (const std::string name : {"MetalBiCGClover", "MetalMGClover"}) {
        std::unique_ptr<LinOpSystemSolver<T>> solver(
            TheLinOpFermSystemSolverFactory::Instance().createObject(
                name, xml, "/InvertParam", state, A));
        T x = zero;
        for (int rhs_index = 0; rhs_index < 2; ++rhs_index) {
          T current = rhs;
          if (rhs_index) {
            source(current, 20261113);
            if (!unpreconditioned) {
              T inverse = zero, hop = zero;
              eo->evenEvenInvLinOp(inverse, current, PLUS);
              eo->oddEvenLinOp(hop, inverse, PLUS);
              current[rb[1]] -= hop;
            }
          }
          // Check a nonzero guess, preserving any outside-subset values.
          source(x, 971 + rhs_index);
          x *= Real(.01);
          auto result = (*solver)(x, current);
          double residual = relative(*A, x, current);
          if (!(residual <= 1e-12))
            throw std::runtime_error("harness FP64 residual failed");
          if (rhs_index == 0 && name == "MetalBiCGClover")
            reference = x;
          double difference = 0;
          if (rhs_index == 0 && name == "MetalMGClover")
            difference = std::sqrt(toDouble(norm2(x - reference, A->subset())) /
                                   toDouble(norm2(reference, A->subset())));
          if (difference > 1e-8)
            throw std::runtime_error("BiCG/MG solution comparison failed");
          if (!unpreconditioned) {
            T even_hop = zero, even_rhs = zero, reconstructed = x;
            eo->evenOddLinOp(even_hop, x, PLUS);
            even_rhs[rb[0]] = b - even_hop;
            if (rhs_index) {
              T source_full = zero;
              source(source_full, 20261113);
              even_rhs[rb[0]] = source_full - even_hop;
            }
            eo->evenEvenInvLinOp(reconstructed, even_rhs, PLUS);
            T source_full = b;
            if (rhs_index)
              source(source_full, 20261113);
            double check = relative(*full, reconstructed, source_full);
            if (check > 1e-12)
              throw std::runtime_error(
                  "full reconstructed Chroma residual failed");
            std::cout << "RECONSTRUCTED solver=" << name << " rhs=" << rhs_index
                      << " full_relative=" << check << '\n';
          }
          std::cout << "HARNESS solver=" << name
                    << " subset=" << (unpreconditioned ? "full" : "odd")
                    << " rhs=" << rhs_index << " K=" << result.n_count
                    << " relative=" << residual
                    << " bicg_mg_difference=" << difference << std::endl;
          // A converged initial guess should require no further iterations.
          auto exact = (*solver)(x, current);
          if (exact.n_count != 0)
            throw std::runtime_error("converged initial guess was not honored");
        }
        T zero_rhs = zero, zero_x = zero;
        auto zero_result = (*solver)(zero_x, zero_rhs);
        if (zero_result.n_count != 0 ||
            toDouble(norm2(zero_x, A->subset())) != 0)
          throw std::runtime_error("zero RHS contract failed");
      }
    }
    // Same array contract as SB::doInversion in the distillation measurement:
    // full spinors -> qprop vector wrapper -> reduced odd Metal vector solve.
    EvenOddPrecCloverFermAct action(create, params);
    std::vector<std::shared_ptr<T>> array_solution(2), array_source(2);
    std::vector<std::shared_ptr<const T>> constant_source(2);
    std::vector<T> array_reference(2);
    int source_time = 176 % dims[3];
    for (int col = 0; col < 2; ++col) {
      array_solution[col] = std::make_shared<T>();
      array_source[col] = std::make_shared<T>();
      source(*array_source[col], 20261114 + col);
      multi1d<int> coord(4);
      size_t volume = 1;
      for (int mu = 0; mu < 4; ++mu)
        volume *= dims[mu];
      for (size_t i = 0; i < volume; ++i) {
        size_t q = i;
        for (int mu = 0; mu < 4; ++mu) {
          coord[mu] = q % dims[mu];
          q /= dims[mu];
        }
        int site = Layout::linearSiteIndex(coord);
        for (int c = 0; c < 3; ++c)
          for (int spin = 0; spin < 4; ++spin)
            if (coord[3] != source_time || spin != 3) {
              auto &z = array_source[col]->elem(site).elem(spin).elem(c);
              z.real() = z.imag() = 0;
            }
      }
      *array_source[col] *= Real(std::sqrt(double(dims[3]) * 4));
      constant_source[col] = array_source[col];
    }
    for (const std::string name : {"MetalBiCGClover", "MetalMGClover"}) {
      GroupXML_t group;
      group.id = name;
      group.xml = inv.str();
      group.path = "/InvertParam";
      std::unique_ptr<SystemSolver<T>> prop(action.qprop(state, group));
      for (auto &solution : array_solution)
        *solution = zero;
      auto results = (*prop)(array_solution, constant_source);
      if (results.size() != 2)
        throw std::runtime_error("qprop array result count");
      for (int col = 0; col < 2; ++col) {
        double residual =
            relative(*full, *array_solution[col], *array_source[col]);
        if (residual > 1e-12)
          throw std::runtime_error("qprop array full FP64 residual failed");
        double difference = 0;
        if (name == "MetalBiCGClover")
          array_reference[col] = *array_solution[col];
        else
          difference = std::sqrt(
              toDouble(norm2(*array_solution[col] - array_reference[col])) /
              toDouble(norm2(array_reference[col])));
        if (difference > 1e-8)
          throw std::runtime_error("qprop array BiCG/MG disagreement");
        std::cout << "QPROP_ARRAY solver=" << name << " rhs=" << col
                  << " source_spin=3 t_source=" << source_time
                  << " full_relative=" << residual
                  << " bicg_mg_difference=" << difference << std::endl;
      }
      std::vector<std::shared_ptr<T>> empty_out;
      std::vector<std::shared_ptr<const T>> empty_in;
      if (!(*prop)(empty_out, empty_in).empty())
        throw std::runtime_error("empty qprop array contract");
    }
#ifdef BUILD_SB
    // Exercise the actual colorvec/spin/tensor adapter used by the user's
    // inline_prop_and_matelem_distillation_superb_w measurement.
    EvenOddPrecCloverFermActEnv::registerAll();
    XMLReader prop_metadata(metadata, "/DBMetaData/Params/Propagator");
    auto action_group = readXMLGroup(prop_metadata, "FermionAction", "FermAct");
    const std::string order = "cxyzXnt";
    auto tensor_dims = SB::latticeSize<Nd + 3>(order, {{'n', 2}, {'t', 1}});
    SB::Tensor<Nd + 3, SB::Complex> colorvec(order, tensor_dims, SB::OnHost,
                                             SB::OnEveryoneReplicated);
    size_t values = 1;
    for (auto length : tensor_dims)
      values *= length;
    double scale = 1 / std::sqrt(double(dims[0]) * dims[1] * dims[2] * 3);
    for (size_t index = 0; index < values; ++index) {
      SB::Coor<Nd + 3> coord{};
      size_t q = index;
      for (size_t d = 0; d < coord.size(); ++d) {
        coord[d] = q % tensor_dims[d];
        q /= tensor_dims[d];
      }
      auto seed = hash(index + 987123);
      const SB::Complex phases[] = {
          {scale, 0}, {-scale, 0}, {0, scale}, {0, -scale}};
      colorvec.set(coord, phases[seed % 4]);
    }
    for (const std::string name : {"MetalBiCGClover", "MetalMGClover"}) {
      GroupXML_t group;
      group.id = name;
      group.xml = inv.str();
      group.path = "/InvertParam";
      SB::ChimeraSolver chimera(action_group, group, u);
      auto solved = SB::doInversion(chimera, colorvec, source_time, 0, dims[3],
                                    {3}, 1, "cxyzXnSst");
      for (int col = 0; col < 2; ++col) {
        T source_field = zero, solution_field = zero;
        colorvec.kvslice_from_size({{'n', col}}, {{'n', 1}})
            .rename_dims({{'n', 's'}})
            .copyTo(SB::asTensorView(source_field)
                        .toComplex()
                        .kvslice_from_size({{'s', 3}, {'t', source_time}},
                                           {{'s', 1}, {'t', 1}}));
        solved.kvslice_from_size({{'n', col}, {'s', 0}}, {{'n', 1}, {'s', 1}})
            .rename_dims({{'s', '!'}, {'S', 's'}})
            .copyTo(SB::asTensorView(solution_field).toComplex());
        double residual = relative(*full, solution_field, source_field);
        if (!(residual <= 1e-12))
          throw std::runtime_error(
              "SB distillation array full FP64 residual failed");
        std::cout << "SUPERB_ARRAY solver=" << name << " rhs=" << col
                  << " max_rhs=1 source_spin=3 t_source=" << source_time
                  << " full_relative=" << residual << std::endl;
      }
    }
#endif
    std::cout
        << "PASS MetalBiCGClover and MetalMGClover: real LIME gauge, native "
           "stout state, "
           "full/odd solves, two reused RHS, nonzero/exact guesses, zero RHS, "
           "native reconstruction, original CPU FP64 residuals <=1e-12"
        << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR " << e.what() << std::endl;
    status = 1;
  } catch (const std::string &e) {
    std::cerr << "ERROR " << e << std::endl;
    status = 1;
  }
  QDP_finalize();
  return status;
}
