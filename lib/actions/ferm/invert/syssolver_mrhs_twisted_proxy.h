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
#include <string>
namespace Chroma {

 namespace MRHSUtils {

	void resetTwist(XMLReader& xml, const Real& Twist)
	{
		std::string sub_solver_type;
		read(xml, "invType", sub_solver_type);

		if( sub_solver_type == "QUDA_WILSON_INVERTER"
				|| sub_solver_type == "QUDA_MULTIGRID_WILSON_INVERTER" ) {
			QDPIO::cout << "QUDA Wilson is not supported in Twisted Linop, only Clover" <<std::endl;
			QDP_abort(1);
		}

		if( sub_solver_type.find("QOP") != std::string::npos) {
			QDPIO::cout << "QOP solver is not supported in Twisted Linop, only Clover" <<std::endl;
			QDP_abort(1);
		}

		if ( sub_solver_type.find("QPHIX") != std::string::npos
				|| sub_solver_type.find("MG_PROTO") != std::string::npos ) {
			QDPIO::cout << "QPHIX and MGPROTO solvers are not supported in Twisted Linop, only Clover" <<std::endl;
			QDP_abort(1);
		}

		if ( sub_solver_type == "QUDA_CLOVER_INVERTER"
				|| sub_solver_type == "QUDA_MULTIGRID_CLOVER_INVERTER"
			    || sub_solver_type.find("MP_CLOVER") != std::string::npos ) {

			const std::string twist_path = "./CloverParams/TwistedM";
			// const std::string rsd_target_path = "./RsdTarget";

			QDPIO::cout << "Resetting Twisted Mass to " << Twist << std::endl;

			if( xml.count(twist_path) != 1) {
				QDPIO::cout << "Error: CloverParams does not contain TwistedMass. Cant reset" <<std::endl;
				QDP_abort(1);
			}

			xml.set<Real>(twist_path, Twist);
			// xml.set<Real>(rsd_target_path, RsdTarget);

			QDPIO::cout << "SUCCESS" << std::endl;
			return;
		}

		// If we haven't returned we are not using
		QDPIO::cout << "Syssolver doesnt need param reset" << std::endl;
		return;
	}
}
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
			const FermAct4D<T,P,Q>& S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_params(params),
				_M_single(static_cast<LinOp<T,P,Q>*>(S_ferm.linOp(state))),
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

			QDPIO::cout << "Setting twist " << _params.Twists[i] << std::endl;
			try {
			MRHSUtils::resetTwist(sub_solver_xml, _params.Twists[i]);
			}
			catch( const std::string& e) {
				QDPIO::cout << "Caught Exception in resetTwist: " << e << std::endl;
				QDP_abort(1);
			}

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
	const Handle<LinOp<T,P,Q>> _M_single;
	const Handle< FermState<T,P,Q> > _state;
};


template<typename T,typename P, typename Q,
	template<typename,typename,typename> class LinOp>
class MdagMMRHSSysSolverTwistedProxy : public MdagMMRHSSystemSolver<T> {
public:

	MdagMMRHSSysSolverTwistedProxy(const SysSolverMRHSTwistedParams& params,
			const FermAct4D<T,P,Q>& S_ferm,
			const Handle< FermState<T,P,Q> > state) :
				_params(params),
				_M_single(static_cast<LinOp<T,P,Q>*>(S_ferm.linOp(state))),
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
			QDPIO::cout << "Setting twist " << _params.Twists[i] << std::endl;
			try {
				MRHSUtils::resetTwist(sub_solver_xml, _params.Twists[i]);
			}
			catch( const std::string& e) {
				QDPIO::cout << "Caught Exception in resetTwist: " << e << std::endl;
				QDP_abort(1);
			}

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
