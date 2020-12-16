/*
 * syssolver_linop_clover_mg_proto.h
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_LINOP_CLOVER_MG_PROTO_QPHIX_H_
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_LINOP_CLOVER_MG_PROTO_QPHIX_H_

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/mg_proto/mgproto_solver_params.h"
#include "actions/ferm/invert/syssolver_linop.h"
using namespace QDP;

namespace Chroma {

//! Registration and other yuckies
  namespace LinOpSysSolverMGProtoQPhiXCloverEnv
  {
    //! Register the syssolver
    bool registerAll();


  }

  class LinOpSysSolverMGProtoQPhiXClover : public LinOpSystemSolver<LatticeFermion>
  {
  public:
	  using T = LatticeFermion;
	  using Q = multi1d<LatticeColorMatrix>;

	  LinOpSysSolverMGProtoQPhiXClover(Handle< LinearOperator<T> > A_,
			  Handle< FermState<T,Q,Q> > state_,
			  const MGProtoSolverParams& invParam_);

	  ~LinOpSysSolverMGProtoQPhiXClover();


	  //! Return the subset on which the operator acts
	  const Subset& subset() const;

	  //! Solve It!
	  SystemSolverResults_t operator()(T& psi, const T& chi) const;

  private:
	  Handle< LinearOperator< T > > A;
	  Handle< FermState<T,Q,Q> > state;
	  MGProtoSolverParams invParam;
  };

};




#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_SYSSOLVER_LINOP_CLOVER_MG_PROTO_H_ */
