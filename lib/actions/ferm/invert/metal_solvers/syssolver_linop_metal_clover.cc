#include "actions/ferm/invert/metal_solvers/syssolver_linop_metal_clover.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace Chroma {
  using namespace MetalCloverEnv;
  namespace {
    struct Complex64 {
      double re, im;
    };
    void check(int status) {
      if (status)
        throw std::runtime_error(metalclover_last_error());
    }
    template <class T>
    void optional(XMLReader &xml, const char *tag, T &value) {
      if (xml.count(tag))
        read(xml, tag, value);
    }
  } // namespace
  SysSolverMetalCloverParams::SysSolverMetalCloverParams(
      XMLReader &xml, const std::string &path)
      : clover(xml, path + "/CloverParams") {
    metalclover_default_options(&options);
    XMLReader top(xml, path);
    Real target;
    read(top, "RsdTarget", target);
    tolerance = toDouble(target);
    optional(top, "MaxIter", max_iterations);
    optional(top, "Workers", options.workers);
    optional(top, "Restart", options.restart);
    optional(top, "NullVectors", options.null_candidates);
    optional(top, "NullIterations", options.null_iterations);
    optional(top, "NullRestart", options.null_restart);
    optional(top, "NullWarm", options.null_warm);
    optional(top, "NullBatch", options.null_batch);
    optional(top, "SecondFilter", options.second_filter);
    optional(top, "BottomIterations", options.bottom_iterations);
    if (top.count("NullTolerance")) {
      Real v;
      read(top, "NullTolerance", v);
      options.null_tolerance = toDouble(v);
    }
    if (top.count("Blocks")) {
      multi1d<int> b;
      read(top, "Blocks", b);
      if (b.size() != 4)
        throw std::runtime_error("Metal Clover requires four Blocks");
      for (int mu = 0; mu < 4; ++mu) {
        if (b[mu] < 1)
          throw std::runtime_error("Metal Clover block must be positive");
        options.blocks[mu] = b[mu];
      }
    }
    if (toDouble(clover.u0) != 1)
      throw std::runtime_error("Metal Clover currently requires u0=1");
    options.mass = toDouble(clover.Mass);
    bool aniso = clover.anisoParam.anisoP;
    if (aniso && clover.anisoParam.t_dir != 3)
      throw std::runtime_error("Metal Clover currently requires t_dir=3");
    if (clover.twisted_m_usedP || clover.max_norm_usedP)
      throw std::runtime_error(
          "Metal Clover does not support TwistedM or MaxNorm");
    options.xi0 = aniso ? toDouble(clover.anisoParam.xi_0) : 1;
    options.nu = aniso ? toDouble(clover.anisoParam.nu) : 1;
    options.clover_r = toDouble(clover.clovCoeffR);
    options.clover_t = aniso ? toDouble(clover.clovCoeffT) : options.clover_r;
    // FermState supplies already smeared, already boundary-modified links.
    options.apply_temporal_antiperiodic = 0;
    if (!(tolerance > 0 && tolerance < 1) || max_iterations < 1)
      throw std::runtime_error("invalid Metal Clover stopping parameters");
  }

  LinOpSysSolverMetalClover::LinOpSysSolverMetalClover(
      Handle<LinearOperator<T>> A_, Handle<FermState<T, Q, Q>> state,
      const SysSolverMetalCloverParams &p, int method)
      : A(A_), params(p),
        name(method == METALCLOVER_MG ? "MetalMGClover" : "MetalBiCGClover") {
#if BASE_PRECISION != 64
    throw std::runtime_error(
        "Metal Clover native residual checks require FP64 Chroma");
#endif
    if (Nd != 4 || Nc != 3 || Ns != 4 || Layout::numNodes() != 1)
      throw std::runtime_error(
          "Metal Clover currently supports single-node 4D Nc=3 Ns=4 QDP++");
    full = A->subset().numSiteTable() == all.numSiteTable();
    if (!full && A->subset().numSiteTable() != rb[1].numSiteTable())
      throw std::runtime_error("Metal Clover requires full or odd subset");
    const auto &dims = Layout::lattSize();
    size_t volume = 1;
    for (int mu = 0; mu < 4; ++mu) {
      params.options.shape[mu] = dims[mu];
      volume *= dims[mu];
    }
    params.options.method = method;
    std::vector<int> lex(volume);
    multi1d<int> coord(4);
    for (size_t i = 0; i < volume; ++i) {
      size_t q = i;
      int parity = 0;
      for (int mu = 0; mu < 4; ++mu) {
        coord[mu] = q % dims[mu];
        q /= dims[mu];
        parity += coord[mu];
      }
      lex[i] = Layout::linearSiteIndex(coord);
      if (full || parity % 2)
        sites.push_back(lex[i]);
    }
    // Require exact subset membership; a generic half-size subset is not
    // enough.
    std::vector<bool> active(volume);
    const auto &table = A->subset().siteTable();
    for (int i = 0; i < table.size(); ++i)
      active[table[i]] = true;
    for (int site : sites)
      if (!active[site])
        throw std::runtime_error("Metal Clover subset is not canonical odd");
    std::vector<Complex64> links(volume * 36);
    const Q &u = state->getLinks();
    for (size_t i = 0; i < volume; ++i)
      for (int mu = 0; mu < 4; ++mu)
        for (int row = 0; row < 3; ++row)
          for (int col = 0; col < 3; ++col) {
            auto z = u[mu].elem(lex[i]).elem().elem(row, col);
            links[(i * 4 + mu) * 9 + row * 3 + col] = {double(z.real()),
                                                       double(z.imag())};
          }
    metalclover_context *created = nullptr;
    check(metalclover_create(&params.options, links.data(), &created));
    std::unique_ptr<metalclover_context, decltype(&metalclover_destroy)> guard(
        created, metalclover_destroy);
    // Compare the supplied Chroma operator with the backend before solving any
    // RHS.
    std::vector<Complex64> probe(sites.size() * 12), native(probe.size()),
        backend(probe.size());
    for (size_t j = 0; j < probe.size(); ++j)
      probe[j] = {std::sin(double(j) * .731 + 1),
                  std::cos(double(j) * .419 + 2)};
    if (apply(this, probe.data(), native.data(), probe.size()))
      throw std::runtime_error("Metal Clover native operator probe failed");
    check(full ? metalclover_apply_full(created, probe.data(), backend.data())
               : metalclover_apply_odd(created, probe.data(), backend.data()));
    double error = 0, reference = 0;
    for (size_t j = 0; j < probe.size(); ++j) {
      double re = native[j].re - backend[j].re,
             im = native[j].im - backend[j].im;
      error += re * re + im * im;
      reference += native[j].re * native[j].re + native[j].im * native[j].im;
    }
    double relative = std::sqrt(error / (reference > 0 ? reference : 1));
    if (!std::isfinite(relative) || relative > 1e-12)
      throw std::runtime_error("Metal Clover differs from supplied Chroma "
                               "operator (action/subset conventions)");
    context = guard.release();
    QDPIO::cout << name << ": subset=" << (full ? "full" : "odd Schur")
                << " operator_match=" << relative
                << " setup_seconds=" << metalclover_setup_seconds(context)
                << std::endl;
  }
  LinOpSysSolverMetalClover::~LinOpSysSolverMetalClover() {
    metalclover_destroy(context);
  }
  void LinOpSysSolverMetalClover::pack(const T &field, void *buffer) const {
    auto *bytes = static_cast<char *>(buffer);
    for (size_t i = 0; i < sites.size(); ++i)
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto z = field.elem(sites[i]).elem(s).elem(c);
          Complex64 v{double(z.real()), double(z.imag())};
          std::memcpy(bytes + (i * 12 + c * 4 + s) * sizeof(v), &v, sizeof(v));
        }
  }
  void LinOpSysSolverMetalClover::unpack(T &field, const void *buffer) const {
    auto *bytes = static_cast<const char *>(buffer);
    for (size_t i = 0; i < sites.size(); ++i)
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          Complex64 v;
          std::memcpy(&v, bytes + (i * 12 + c * 4 + s) * sizeof(v), sizeof(v));
          auto &z = field.elem(sites[i]).elem(s).elem(c);
          z.real() = v.re;
          z.imag() = v.im;
        }
  }
  int LinOpSysSolverMetalClover::apply(void *self, const void *in, void *out,
                                       size_t n) {
    try {
      auto &solver = *static_cast<LinOpSysSolverMetalClover *>(self);
      if (n != solver.sites.size() * 12)
        return 1;
      T x = zero, y = zero;
      solver.unpack(x, in);
      (*solver.A)(y, x, PLUS);
      solver.pack(y, out);
      return 0;
    } catch (...) {
      return 1;
    }
  }
  SystemSolverResults_t
  LinOpSysSolverMetalClover::operator()(T &psi, const T &chi) const {
    try {
      StopWatch total;
      total.start();
      std::vector<Complex64> b(sites.size() * 12), x(b.size());
      pack(chi, b.data());
      pack(psi, x.data());
      auto *self = const_cast<LinOpSysSolverMetalClover *>(this);
      int status =
          full ? metalclover_solve_full(context, b.data(), x.data(),
                                        params.tolerance, params.max_iterations,
                                        apply, self, &last_info)
               : metalclover_solve_odd(context, b.data(), x.data(),
                                       params.tolerance, params.max_iterations,
                                       apply, self, &last_info);
      check(status);
      unpack(psi, x.data());
      StopWatch check_timer;
      check_timer.start();
      T residual = zero;
      (*A)(residual, psi, PLUS);
      residual[subset()] = chi - residual;
      double absolute = std::sqrt(toDouble(norm2(residual, subset())));
      double source = std::sqrt(toDouble(norm2(chi, subset())));
      if (!std::isfinite(absolute) ||
          absolute > params.tolerance * (source > 0 ? source : 1))
        throw std::runtime_error(
            "final original Chroma FP64 residual check failed");
      check_timer.stop();
      total.stop();
      last_info.residual_seconds += check_timer.getTimeInSeconds();
      last_info.seconds = total.getTimeInSeconds();
      last_info.nominal_operator_gflops =
          double(last_info.nominal_operator_flops) / (last_info.seconds * 1e9);
      SystemSolverResults_t result;
      result.n_count = last_info.krylov_iterations;
      result.resid = absolute;
      QDPIO::cout << name << ": iterations=" << result.n_count
                  << " outer=" << last_info.outer_iterations
                  << " relative_residual="
                  << absolute / (source > 0 ? source : 1)
                  << " seconds=" << last_info.seconds
                  << " nominal_operator_GFLOP_per_s="
                  << last_info.nominal_operator_gflops << std::endl;
      return result;
    } catch (const std::exception &e) {
      QDPIO::cerr << name << ": " << e.what() << std::endl;
      QDP_abort(1);
      throw;
    }
  }
  std::vector<SystemSolverResults_t> LinOpSysSolverMetalClover::operator()(
      const std::vector<std::shared_ptr<T>> &solutions,
      const std::vector<std::shared_ptr<const T>> &sources) const {
    try {
      if (solutions.size() != sources.size())
        throw std::runtime_error(
            "Metal Clover source/solution array lengths differ");
      for (size_t i = 0; i < solutions.size(); ++i) {
        if (!solutions[i] || !sources[i])
          throw std::runtime_error("Metal Clover null field in RHS array");
        for (size_t j = 0; j < solutions.size(); ++j)
          if (solutions[i].get() == sources[j].get() ||
              (i != j && solutions[i].get() == solutions[j].get()))
            throw std::runtime_error(
                "Metal Clover array fields must be disjoint");
      }
      last_array_info.clear();
      last_array_info.reserve(sources.size());
      std::vector<SystemSolverResults_t> results;
      results.reserve(sources.size());
      StopWatch watch;
      watch.start();
      uint64_t flops = 0;
      // The caller owns the array; only this column's extra FP64 fields and
      // one GPU workspace are resident in the backend. Hierarchy is reused.
      for (size_t i = 0; i < sources.size(); ++i) {
        results.push_back((*this)(*solutions[i], *sources[i]));
        last_array_info.push_back(last_info);
        flops += last_info.nominal_operator_flops;
      }
      watch.stop();
      QDPIO::cout << name << " array: rhs=" << sources.size()
                  << " resident_GPU_rhs=1 hierarchy_reused=true seconds="
                  << watch.getTimeInSeconds()
                  << " nominal_operator_GFLOP_per_s="
                  << (watch.getTimeInSeconds() > 0
                          ? double(flops) / (1e9 * watch.getTimeInSeconds())
                          : 0)
                  << std::endl;
      return results;
    } catch (const std::exception &e) {
      QDPIO::cerr << name << ": " << e.what() << std::endl;
      QDP_abort(1);
      throw;
    }
  }

  namespace LinOpSysSolverMetalCloverEnv {
    namespace {
      using FS = FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
                           multi1d<LatticeColorMatrix>>;
      LinOpSystemSolver<LatticeFermion> *
      createBiCG(XMLReader &xml, const std::string &path, Handle<FS> state,
                 Handle<LinearOperator<LatticeFermion>> A) {
        return new LinOpSysSolverMetalClover(
            A, state, SysSolverMetalCloverParams(xml, path),
            METALCLOVER_BICGSTAB);
      }
      LinOpSystemSolver<LatticeFermion> *
      createMG(XMLReader &xml, const std::string &path, Handle<FS> state,
               Handle<LinearOperator<LatticeFermion>> A) {
        return new LinOpSysSolverMetalClover(
            A, state, SysSolverMetalCloverParams(xml, path), METALCLOVER_MG);
      }
    } // namespace
    bool registerAll() {
      static bool registered = false;
      bool success = true;
      if (!registered) {
        success &= TheLinOpFermSystemSolverFactory::Instance().registerObject(
            "MetalBiCGClover", createBiCG);
        success &= TheLinOpFermSystemSolverFactory::Instance().registerObject(
            "MetalMGClover", createMG);
        registered = success;
      }
      return success;
    }
  } // namespace LinOpSysSolverMetalCloverEnv
} // namespace Chroma
