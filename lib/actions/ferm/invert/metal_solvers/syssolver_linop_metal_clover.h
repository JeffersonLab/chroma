#ifndef CHROMA_SYSSOLVER_LINOP_METAL_CLOVER_H
#define CHROMA_SYSSOLVER_LINOP_METAL_CLOVER_H
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "handle.h"
#include "state.h"
#include <metalclover/metalclover.h>
#include <vector>

namespace Chroma {
  namespace LinOpSysSolverMetalCloverEnv {
    bool registerAll();
  }

  struct SysSolverMetalCloverParams {
    SysSolverMetalCloverParams(XMLReader &xml, const std::string &path);
    CloverFermActParams clover;
    MetalCloverEnv::metalclover_options options;
    double tolerance = 1e-12;
    int max_iterations = 65536;
  };

  // One serial Metal context (and MG hierarchy) is reused for every call.
  // Supports the unpreconditioned operator and Chroma's asymmetric odd Schur
  // operator A_oo-D_oe A_ee^-1 D_eo. The native PLUS operator checks every
  // refinement residual in CPU FP64.
  class LinOpSysSolverMetalClover : public LinOpSystemSolver<LatticeFermion> {
  public:
    using T = LatticeFermion;
    using Q = multi1d<LatticeColorMatrix>;
    using LinOpSystemSolver<T>::operator();
    LinOpSysSolverMetalClover(Handle<LinearOperator<T>> A,
                              Handle<FermState<T, Q, Q>> state,
                              const SysSolverMetalCloverParams &params,
                              int method);
    ~LinOpSysSolverMetalClover() override;
    const Subset &subset() const override { return A->subset(); }
    SystemSolverResults_t operator()(T &psi, const T &chi) const override;
    std::vector<SystemSolverResults_t> operator()(
        const std::vector<std::shared_ptr<T>> &solutions,
        const std::vector<std::shared_ptr<const T>> &sources) const override;
    const std::vector<MetalCloverEnv::metalclover_solve_info> &
    lastArraySolveInfo() const {
      return last_array_info;
    }
    const MetalCloverEnv::metalclover_solve_info &lastSolveInfo() const {
      return last_info;
    }

  private:
    Handle<LinearOperator<T>> A;
    SysSolverMetalCloverParams params;
    MetalCloverEnv::metalclover_context *context = nullptr;
    bool full = false;
    std::string name;
    std::vector<int> sites;
    mutable MetalCloverEnv::metalclover_solve_info last_info{};
    mutable std::vector<MetalCloverEnv::metalclover_solve_info> last_array_info;
    static int apply(void *self, const void *in, void *out, size_t n);
    void pack(const T &field, void *buffer) const;
    void unpack(T &field, const void *buffer) const;
    LinOpSysSolverMetalClover(const LinOpSysSolverMetalClover &) = delete;
    LinOpSysSolverMetalClover &
    operator=(const LinOpSysSolverMetalClover &) = delete;
  };
} // namespace Chroma
#endif
