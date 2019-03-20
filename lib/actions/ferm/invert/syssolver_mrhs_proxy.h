/*
 * syssolver_mrhs_proxy.h
 *
 *  Created on: Mar 13, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_H_

#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
#include "syssolver.h"
#include "fermact.h"
#include "syssolver_mrhs_proxy_params.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include <sstream>

namespace Chroma {

 //! CG1 system solver namespace
 namespace LinOpSysSolverMRHSProxyEnv
 {
 	 bool registerAll();
 }

 //! CG1 system solver namespace
  namespace MdagMSysSolverMRHSProxyEnv
  {
  	 bool registerAll();
  }

template<typename T,typename P, typename Q>
class LinOpMRHSSysSolverProxy : public LinOpMRHSSystemSolver<T> {
public:
	LinOpMRHSSysSolverProxy(const SysSolverMRHSProxyParams& params,
		    const FermAct4D<T,P,Q>& S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_N(params.BlockSize),
				_M_single(S_ferm.linOp(state)),
				_solver_single(S_ferm.invLinOp(state,params.SubInverterXML))
	{}


	SystemSolverResultsMRHS_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const {

		SystemSolverResultsMRHS_t res(_N);

		// Do the do.
		for(int i=0; i < _N; ++i) {
			SystemSolverResults_t res_tmp = (*_solver_single)(psi[i],chi[i]);
			res.n_count[i] = res_tmp.n_count;
			res.resid[i] = res_tmp.resid;
		}
		return res;
	}

	const Subset& subset() const {
		return _M_single->subset();
	}

	int size() const {
		return _N;
	}
private:
	const int _N;
	const Handle<LinearOperator<T>> _M_single;
	const Handle<SystemSolver<T>> _solver_single;
};

template<typename T,typename P, typename Q>
class MdagMMRHSSysSolverProxy : public MdagMMRHSSystemSolver<T> {
public:
	MdagMMRHSSysSolverProxy(const SysSolverMRHSProxyParams& params,
			const FermAct4D<T,P,Q>& S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_N(params.BlockSize),
				_M_single(S_ferm.linOp(state)),
				_solver_single(S_ferm.invMdagM(state,params.SubInverterXML))
	{}



	SystemSolverResultsMRHS_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const {

		SystemSolverResultsMRHS_t res(_N);

		// Do the do.
		for(int i=0; i < _N; ++i) {
			SystemSolverResults_t res_tmp = (*_solver_single)(psi[i],chi[i]);
			res.n_count[i] = res_tmp.n_count;
			res.resid[i] = res_tmp.resid;
		}
		return res;
	}

	const Subset& subset() const {
		return _M_single->subset();
	}

	int size() const {
		return _N;
	}
private:
	const int _N;
	const Handle<LinearOperator<T>> _M_single;
	Handle<SystemSolver<T>> _solver_single;
};


}

#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_H_ */
