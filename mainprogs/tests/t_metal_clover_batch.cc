// Real-field qprop throughput with one solver/hierarchy reused across RHS
// tiles.
#include "actions/ferm/fermacts/eoprec_clover_fermact_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "util/gauge/reunit.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Chroma;
namespace {
  using T = LatticeFermion;
  using Q = multi1d<LatticeColorMatrix>;
  struct BenchmarkComplex {
    double re, im;
  };
  uint64_t hash(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
  }
  void source(T &x, size_t rhs, const std::vector<int> &sites,
              int source_time) {
    x = zero;
    auto dims = Layout::lattSize();
    size_t spatial = size_t(dims[0]) * dims[1] * dims[2];
    size_t vector = rhs / Ns;
    int source_spin = rhs % Ns;
    double scale = 1 / std::sqrt(double(spatial) * 3);
    for (size_t q = size_t(source_time) * spatial;
         q < size_t(source_time + 1) * spatial; ++q)
      for (int c = 0; c < 3; ++c) {
        auto phase =
            hash((q % spatial) * 3 + c + (20261113ull + vector) * 771337) & 3;
        auto &z = x.elem(sites[q]).elem(source_spin).elem(c);
        z.real() = phase == 0 ? scale : phase == 1 ? -scale : 0;
        z.imag() = phase == 2 ? scale : phase == 3 ? -scale : 0;
      }
  }
  void inverter(XMLBufferWriter &xml, const CloverFermActParams &p,
                const std::string &name) {
    push(xml, "InvertParam");
    write(xml, "invType", name);
    write(xml, "RsdTarget", Real(1e-12));
    write(xml, "MaxIter", 65536);
    write(xml, "RHSBatch", 16);
    multi1d<int> blocks(4);
    blocks[0] = blocks[1] = blocks[2] = 3;
    blocks[3] = 8;
    write(xml, "Blocks", blocks);
    write(xml, "SecondFilter", 16);
    write(xml, "BottomIterations", 4);
    push(xml, "CloverParams");
    write(xml, "Mass", p.Mass);
    write(xml, "u0", p.u0);
    write(xml, "AnisoParam", p.anisoParam);
    write(xml, "clovCoeffR", p.clovCoeffR);
    write(xml, "clovCoeffT", p.clovCoeffT);
    pop(xml);
    pop(xml);
  }
  void extract(const T &x, const std::vector<int> &sites,
               std::vector<BenchmarkComplex> &out) {
    for (size_t i = 0; i < sites.size(); ++i)
      for (int c = 0; c < 3; ++c)
        for (int spin = 0; spin < 4; ++spin) {
          auto z = x.elem(sites[i]).elem(spin).elem(c);
          out[i * 12 + c * 4 + spin] = {double(z.real()), double(z.imag())};
        }
  }
} // namespace
int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: t_metal_clover_batch GAUGE_LIME ACTION_XML "
                 "OUTPUT_PREFIX BASE_NUM_VECS\n";
    return 2;
  }
  std::string gauge = argv[1], action_xml = argv[2], prefix = argv[3];
  int base_vectors = std::stoi(argv[4]);
  if (base_vectors < 1 || base_vectors > 24)
    throw std::runtime_error(
        "BASE_NUM_VECS must be 1..24 (4..96 physical RHS)");
  int base = Ns * base_vectors;
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    static_assert(BASE_PRECISION == 64, "requires native FP64 Chroma");
    XMLReader metadata(action_xml);
    multi1d<int> dims;
    read(metadata, "/DBMetaData/lattSize", dims);
    if (dims.size() != 4 || dims[0] != 12 || dims[1] != 12 || dims[2] != 12 ||
        dims[3] != 96)
      throw std::runtime_error("this performance harness targets 12^3x96");
    Layout::setLattSize(dims);
    Layout::create();
    if (Layout::numNodes() != 1)
      throw std::runtime_error("one MPI rank required");
    std::vector<int> sites(size_t(12) * 12 * 12 * 96);
    multi1d<int> coord(4);
    for (size_t i = 0; i < sites.size(); ++i) {
      size_t q = i;
      for (int mu = 0; mu < 4; ++mu) {
        coord[mu] = q % dims[mu];
        q /= dims[mu];
      }
      sites[i] = Layout::linearSiteIndex(coord);
    }
    const std::string path = "/DBMetaData/Params/Propagator/FermionAction";
    CloverFermActParams params(metadata, path);
    Q u(4);
    XMLReader file_xml, record_xml;
    QDPFileReader reader(file_xml, gauge, QDPIO_SERIAL);
    read(reader, record_xml, u);
    if (reader.bad())
      throw std::runtime_error("gauge LIME/checksum read failed");
    reader.close();
    for (int mu = 0; mu < 4; ++mu)
      reunit(u[mu]);
    CreateFermStateEnv::registerAll();
    LinOpSysSolverEnv::registerAll();
    auto create = CreateFermStateEnv::reader(metadata, path);
    Handle<FermState<T, Q, Q>> state((*create)(u));
    UnprecCloverLinOp full(state, params);
    EvenOddPrecCloverFermAct action(create, params);
    std::string reference_path = prefix + ".reference.bin";
    if (std::ifstream(reference_path).good() ||
        std::ifstream(prefix + ".csv").good())
      throw std::runtime_error(
          "refusing to replace existing benchmark outputs");
    std::ofstream csv(prefix + ".csv");
    if (!csv)
      throw std::runtime_error("open output CSV");
    csv << std::setprecision(17)
        << "method,num_vecs,source_spins,rhs,caller_tile,seconds,seconds_per_"
           "rhs,"
           "seconds_per_colorvec,rhs_per_second,max_"
           "full_relative,max_solution_difference,min_iterations,max_"
           "iterations,total_iterations\n";
    std::cout << std::setprecision(17);
    const int source_time = 176 % dims[3];
    std::cout << "DISTILLATION_SOURCES num_vecs=" << base_vectors
              << " source_spins=" << Ns << " physical_rhs=" << base
              << " t_source=" << source_time
              << " seed_base=20261113 spin_fast=true" << std::endl;
    bool have_reference = false;
    for (const std::string name : {"MetalBiCGClover", "MetalMGClover"}) {
      XMLBufferWriter inv;
      inverter(inv, params, name);
      GroupXML_t group;
      group.id = name;
      group.xml = inv.str();
      group.path = "/InvertParam";
      StopWatch setup;
      setup.start();
      Handle<SystemSolver<T>> prop(action.qprop(state, group));
      setup.stop();
      std::cout << "BATCH_SETUP method=" << name
                << " seconds=" << setup.getTimeInSeconds()
                << " hierarchy_instances=" << (name == "MetalMGClover" ? 1 : 0)
                << std::endl;
      std::vector<std::shared_ptr<T>> solutions;
      std::vector<std::shared_ptr<const T>> sources;
      auto extend = [&](int count) {
        while (int(sources.size()) < count) {
          auto b = std::make_shared<T>();
          source(*b, sources.size(), sites, source_time);
          sources.push_back(b);
          solutions.push_back(std::make_shared<T>(zero));
        }
      };
      for (auto test : {std::pair<int, int>{base, 1},
                        {base, 8},
                        {base, 16},
                        {2 * base, 16}}) {
        int count = test.first, tile = test.second;
        extend(count);
        for (int i = 0; i < count; ++i)
          *solutions[i] = zero;
        std::vector<SystemSolverResults_t> results(count);
        std::cout << "BATCH_BEGIN method=" << name << " num_vecs=" << count / Ns
                  << " source_spins=" << Ns << " rhs=" << count
                  << " caller_tile=" << tile << std::endl;
        StopWatch watch;
        watch.start();
        for (int i = 0; i < count; i += tile) {
          int end = std::min(count, i + tile);
          auto out = std::vector<std::shared_ptr<T>>(solutions.begin() + i,
                                                     solutions.begin() + end);
          auto in = std::vector<std::shared_ptr<const T>>(
              sources.begin() + i, sources.begin() + end);
          auto result = (*prop)(out, in);
          if (result.size() != size_t(end - i))
            throw std::runtime_error("array result count mismatch");
          std::copy(result.begin(), result.end(), results.begin() + i);
        }
        watch.stop();
        double max_residual = 0, max_difference = 0;
        int min_iterations = 65536, max_iterations = 0;
        uint64_t iterations = 0;
        std::vector<BenchmarkComplex> field(sites.size() * 12),
            reference(field.size());
        std::ofstream save;
        std::ifstream compare;
        if (!have_reference)
          save.open(reference_path, std::ios::binary);
        else
          compare.open(reference_path, std::ios::binary);
        if ((!have_reference && !save) || (have_reference && !compare))
          throw std::runtime_error("open scalar reference fields");
        for (int i = 0; i < count; ++i) {
          T r = zero;
          full(r, *solutions[i], PLUS);
          r -= *sources[i];
          double bnorm2 = toDouble(norm2(*sources[i]));
          double relative = std::sqrt(toDouble(norm2(r)) / bnorm2);
          double wrapper_relative =
              toDouble(results[i].resid) / std::sqrt(bnorm2);
          if (!(relative <= 1e-12) || !(wrapper_relative <= 1e-12))
            throw std::runtime_error(
                "full original Chroma FP64 residual failed");
          max_residual = std::max(max_residual, relative);
          min_iterations = std::min(min_iterations, results[i].n_count);
          max_iterations = std::max(max_iterations, results[i].n_count);
          iterations += results[i].n_count;
          if (i < base) {
            extract(*solutions[i], sites, field);
            if (!have_reference) {
              save.write(reinterpret_cast<const char *>(field.data()),
                         field.size() * sizeof(BenchmarkComplex));
              if (!save)
                throw std::runtime_error("write scalar reference");
            } else {
              compare.read(reinterpret_cast<char *>(reference.data()),
                           reference.size() * sizeof(BenchmarkComplex));
              if (!compare)
                throw std::runtime_error("read scalar reference");
              double error = 0, norm = 0;
              for (size_t q = 0; q < field.size(); ++q) {
                double re = field[q].re - reference[q].re,
                       im = field[q].im - reference[q].im;
                error += re * re + im * im;
                norm += reference[q].re * reference[q].re +
                        reference[q].im * reference[q].im;
              }
              double difference = std::sqrt(error / norm);
              if (!(difference <= 1e-8))
                throw std::runtime_error(
                    "scalar/batch/BiCG/MG solution disagreement");
              max_difference = std::max(max_difference, difference);
            }
          }
        }
        have_reference = true;
        double elapsed = watch.getTimeInSeconds();
        csv << name << ',' << count / Ns << ',' << Ns << ',' << count << ','
            << tile << ',' << elapsed << ',' << elapsed / count << ','
            << elapsed / (count / Ns) << ',' << count / elapsed << ','
            << max_residual << ',' << max_difference << ',' << min_iterations
            << ',' << max_iterations << ',' << iterations << '\n';
        csv.flush();
        std::cout << "BATCH_RESULT method=" << name
                  << " num_vecs=" << count / Ns << " source_spins=" << Ns
                  << " rhs=" << count << " caller_tile=" << tile
                  << " seconds=" << elapsed
                  << " seconds_per_rhs=" << elapsed / count
                  << " seconds_per_colorvec=" << elapsed / (count / Ns)
                  << " max_full_relative=" << max_residual
                  << " max_solution_difference=" << max_difference
                  << " K_min=" << min_iterations << " K_max=" << max_iterations
                  << std::endl;
      }
    }
    std::cout << "PASS real 12^3x96 four-spin distillation qprop "
                 "scalar/batched RHS, one reused "
                 "hierarchy, all full FP64 residuals <=1e-12"
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
