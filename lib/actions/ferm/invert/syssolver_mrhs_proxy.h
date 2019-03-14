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
template<typename T,typename P, typename Q>
class LinOpMRHSSysSolverProxy : public LinOpMRHSSystemSolver<T> {
public:
	LinOpMRHSSysSolverProxy(const SysSolverMRHSProxyParams& params,
			const Handle< FermAct4D<T,P,Q> > S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_params(params),
				_M_single(S_ferm->linOp(state))
	{
		std::istringstream sub_solver_xml_stream(_params.SubInverterXML.xml);
		XMLReader sub_solver_xml(sub_solver_xml_stream);

		QDPIO::cout << "Creating SubSolver " << std::endl;
		_solver_single = TheLinOpFermSystemSolverFactory::Instance().createObject(
				_params.SubInverterXML.id,
				sub_solver_xml,
				_params.SubInverterXML.path,
				state,
				_M_single);
		QDPIO::cout << "Done" << std::endl;
	}


	SystemSolverResultsMRHS_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const {
		int N = size();
		SystemSolverResultsMRHS_t res(N);

		// Do the do.
		for(int i=0; i < N; ++i) {
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
		return _params.BlockSize;
	}
private:
	const SysSolverMRHSProxyParams _params;
	const Handle<LinearOperator<T>> _M_single;
	Handle<SystemSolver<T>> _solver_single;
};

template<typename T,typename P, typename Q>
class MdagMMRHSSysSolverProxy : public MdagMMRHSSystemSolver<T> {
public:
	MdagMMRHSSysSolverProxy(const SysSolverMRHSProxyParams& params,
			const Handle< FermAct4D<T,P,Q> > S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_params(params),
				_M_single(S_ferm->linOp(state))
	{
		std::istringstream sub_solver_xml_stream(_params.SubInverterXML.xml);
		XMLReader sub_solver_xml(sub_solver_xml_stream);

		QDPIO::cout << "Creating SubSolver " << std::endl;
		_solver_single = TheMdagMFermSystemSolverFactory::Instance().createObject(
				_params.SubInverterXML.id,
				sub_solver_xml,
				_params.SubInverterXML.path,
				state,
				_M_single);
		QDPIO::cout << "Done" << std::endl;
	}


	SystemSolverResultsMRHS_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const {
		int N = size();
		SystemSolverResultsMRHS_t res(N);

		// Do the do.
		for(int i=0; i < N; ++i) {
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
		return _params.BlockSize;
	}
private:
	const SysSolverMRHSProxyParams _params;
	const Handle<LinearOperator<T>> _M_single;
	Handle<SystemSolver<T>> _solver_single;
};


}

#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_PROXY_H_ */
