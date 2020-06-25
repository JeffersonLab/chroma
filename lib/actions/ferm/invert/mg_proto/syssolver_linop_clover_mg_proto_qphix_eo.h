/*
 * syssolver_linop_clover_mg_proto.h
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_LINOP_CLOVER_MG_PROTO_QPHIX_EO_H_
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_LINOP_CLOVER_MG_PROTO_QPHIX_EO_H_

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
  namespace LinOpSysSolverMGProtoQPhiXEOCloverEnv
  {
    //! Register the syssolver
    bool registerAll();


  }

  using EoFGMRES = const MG::FGMRESSolverQPhiX;

  class LinOpSysSolverMGProtoQPhiXEOClover : public LinOpSystemSolver<LatticeFermion>
  {
  public:
	  using T = LatticeFermion;
	  using Q = multi1d<LatticeColorMatrix>;

	  LinOpSysSolverMGProtoQPhiXEOClover(Handle< LinearOperator<T> > A_,
			  Handle< FermState<T,Q,Q> > state_,
			  const MGProtoSolverParams& invParam_);

	  ~LinOpSysSolverMGProtoQPhiXEOClover();


	  //! Return the subset on which the operator acts
	  const Subset& subset() const;

	  //! Solve It!
	  SystemSolverResults_t operator()(T& psi, const T& chi) const override;
  	  std::vector<SystemSolverResults_t> operator()(const std::vector<std::shared_ptr<T>>& psi, const std::vector<std::shared_ptr<const T>>& chi) const override;

  private:
	  Handle< LinearOperator< T > > A;
	  Handle< FermState<T,Q,Q> > state;
	  MGProtoSolverParams invParam;
	  const std::string subspaceId;
	  std::shared_ptr<MGProtoHelpersQPhiX::MGPreconditionerEO> mg_pointer;
	  std::shared_ptr<MG::QPhiXWilsonCloverEOLinearOperator > M_ptr;

	  // Shorthand for the UnprecWrapper
	  using UnprecFGMRES =  MG::UnprecFGMRESSolverQPhiXWrapper;


	  std::shared_ptr<UnprecFGMRES> wrapped;
	  std::shared_ptr<EoFGMRES> eo_solver;

  };

};




#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_LINOP_CLOVER_MG_PROTO_H_ */
