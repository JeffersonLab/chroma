/*
 * Solve MdagM*psi=chi system by mg_proto
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_MDAGM_CLOVER_MG_PROTO_QPHIX_EO_H_
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_MDAGM_CLOVER_MG_PROTO_QPHIX_EO_H_

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/mg_proto/mgproto_solver_params.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/mg_proto/mg_proto_qphix_helpers.h"
#include "lattice/solver.h"
#include "lattice/fgmres_common.h"
#include "lattice/qphix/invfgmres_qphix.h"
#include "lattice/qphix/qphix_qdp_utils.h"
#include "lattice/qphix/qphix_eo_clover_linear_operator.h"
#include "actions/ferm/invert/mg_solver_exception.h"
#include <memory>

using namespace QDP;

namespace Chroma {

//! Registration and other yuckies
  namespace MdagMSysSolverMGProtoQPhiXEOCloverEnv
  {
    //! Register the syssolver
    bool registerAll();


  }

  using EoFGMRES = const MG::FGMRESSolverQPhiX;

  class MdagMSysSolverMGProtoQPhiXEOClover : public MdagMSystemSolver<LatticeFermion>
  {
  public:
	  using T = LatticeFermion;
	  using Q = multi1d<LatticeColorMatrix>;

	  MdagMSysSolverMGProtoQPhiXEOClover(Handle< LinearOperator<T> > A_,
			  Handle< FermState<T,Q,Q> > state_,
			  const MGProtoSolverParams& invParam_);

	  ~MdagMSysSolverMGProtoQPhiXEOClover();


	  //! Return the subset on which the operator acts
	  const Subset& subset() const;

	  //! Solve It!
	  SystemSolverResults_t operator()(T& psi, const T& chi) const override;

	  SystemSolverResults_t operator()(T& psi, const T& chi,
              AbsChronologicalPredictor4D<T>& predictor) const override;
  private:
	  Handle< LinearOperator< T > > A;
	  Handle< FermState<T,Q,Q> > state;
	  MGProtoSolverParams invParam;
      MG::FGMRESParams fine_solve_params;
      std::string subspaceId;
	  mutable std::shared_ptr<MGProtoHelpersQPhiX::MGPreconditionerEO> mg_pointer;
	  mutable std::shared_ptr<MG::QPhiXWilsonCloverEOLinearOperator > M_ptr;

	  // Shorthand for the UnprecWrapper
	  using UnprecFGMRES =  MG::UnprecFGMRESSolverQPhiXWrapper;
	  mutable std::shared_ptr<EoFGMRES> eo_solver;

  };

};




#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_MDAGM_CLOVER_MG_PROTO_H_ */
