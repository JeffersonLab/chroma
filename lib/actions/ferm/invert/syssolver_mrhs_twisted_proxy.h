/*
 * syssolver_mrhs_twisted_proxy.h
 *
 *  Created on: Mar 15, 2019
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_TWISTED_PROXY_H_
#define LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_TWISTED_PROXY_H_

#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
#include "syssolver.h"
#include "fermact.h"
#include "syssolver_mrhs_twisted_params.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/linop/shifted_linop_w.h"
#include <sstream>

namespace Chroma {

 //! CG1 system solver namespace
 namespace LinOpSysSolverMRHSTwistedProxyEnv
 {
 	 bool registerAll();
 }

 //! CG1 system solver namespace
  namespace MdagMSysSolverMRHSTwistedProxyEnv
  {
  	 bool registerAll();
  }

template<typename T,typename P, typename Q,
	template<typename,typename,typename> class LinOp>
class LinOpMRHSSysSolverTwistedProxy : public LinOpMRHSSystemSolver<T> {
public:

	LinOpMRHSSysSolverTwistedProxy(const SysSolverMRHSTwistedParams& params,
			const Handle< FermAct4D<T,P,Q> > S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_params(params),
				_S_ferm(S_ferm),
				_M_single(static_cast<LinOp<T,P,Q>*>(S_ferm->linOp(state))),
				_state(state)
	{}


	SystemSolverResultsMRHS_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const {
		int N = size();
		SystemSolverResultsMRHS_t res(N);

		// Do the do.
		for(int i=0; i < N; ++i) {
			// Make the twisted operator...
			using ShiftedOpT = TwistedShiftedLinOp<T,P,Q,LinOp>;
			Handle<LinearOperator<T>> M_shifted(static_cast<LinearOperator<T>*>(new ShiftedOpT((*_M_single), _params.Twists[i])));

			std::istringstream sub_solver_xml_stream(_params.SubInverterXML.xml);
			XMLReader sub_solver_xml(sub_solver_xml_stream);





			QDPIO::cout << "Creating SubSolver " << std::endl;
			Handle<SystemSolver<T>> _solver_single = TheLinOpFermSystemSolverFactory::Instance().createObject(
					_params.SubInverterXML.id,
					sub_solver_xml,
					_params.SubInverterXML.path,
					_state,
					M_shifted);


			QDPIO::cout << "Done" << std::endl;

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
	const SysSolverMRHSTwistedParams _params;
	const Handle< FermAct4D<T,P,Q> > _S_ferm;
	const Handle<LinOp<T,P,Q>> _M_single;
	const Handle< FermState<T,P,Q> > _state;
};


template<typename T,typename P, typename Q,
	template<typename,typename,typename> class LinOp>
class MdagMMRHSSysSolverTwistedProxy : public MdagMMRHSSystemSolver<T> {
public:

	MdagMMRHSSysSolverTwistedProxy(const SysSolverMRHSTwistedParams& params,
			const Handle< FermAct4D<T,P,Q> > S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_params(params),
				_S_ferm(S_ferm),
				_M_single(static_cast<LinOp<T,P,Q>*>(S_ferm->linOp(state))),
				_state(state)
	{}


	SystemSolverResultsMRHS_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const {
		int N = size();
		SystemSolverResultsMRHS_t res(N);

		// Do the do.
		for(int i=0; i < N; ++i) {
			// Make the twisted operator...
			using ShiftedOpT = TwistedShiftedLinOp<T,P,Q,LinOp>;
			Handle<LinearOperator<T>> M_shifted(static_cast<LinearOperator<T>*>(new ShiftedOpT((*_M_single), _params.Twists[i])));

			std::istringstream sub_solver_xml_stream(_params.SubInverterXML.xml);
			XMLReader sub_solver_xml(sub_solver_xml_stream);


			QDPIO::cout << "Creating SubSolver " << std::endl;
			Handle<SystemSolver<T>> _solver_single = TheMdagMFermSystemSolverFactory::Instance().createObject(
					_params.SubInverterXML.id,
					sub_solver_xml,
					_params.SubInverterXML.path,
					_state,
					M_shifted);


			QDPIO::cout << "Done" << std::endl;

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
	const SysSolverMRHSTwistedParams _params;
	const Handle< FermAct4D<T,P,Q> > _S_ferm;
	const Handle<LinOp<T,P,Q>> _M_single;
	const Handle< FermState<T,P,Q> > _state;
};

template<typename T, typename P, typename Q>
using SymEvenOddPrecLogDetLinOpMRHSSysSolverTwistedProxy =
		LinOpMRHSSysSolverTwistedProxy<T,P,Q,SymEvenOddPrecLogDetLinearOperator>;

template<typename T, typename P, typename Q>
using EvenOddPrecLinOpMRHSSysSolverTwistedProxy =
		LinOpMRHSSysSolverTwistedProxy<T,P,Q,EvenOddPrecLinearOperator>;



template<typename T, typename P, typename Q>
using SymEvenOddPrecLogDetMdagMMRHSSysSolverTwistedProxy =
		MdagMMRHSSysSolverTwistedProxy<T,P,Q,SymEvenOddPrecLogDetLinearOperator>;

template<typename T, typename P, typename Q>
using EvenOddPrecMdagMMRHSSysSolverTwistedProxy =
		MdagMMRHSSysSolverTwistedProxy<T,P,Q,EvenOddPrecLinearOperator>;

}



#endif /* LIB_ACTIONS_FERM_INVERT_SYSSOLVER_MRHS_TWISTED_PROXY_H_ */
