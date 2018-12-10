/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "update/molecdyn/monomial/seoprec_constdet_two_flavor_multihasen_cancle_monomial_w.h"

namespace Chroma
{
	namespace SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomialEnv
	{
		namespace
		{
			//! Callback function for the factory
			Monomial<multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const std::string& path)
				{
					return new SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(
							TwoFlavorMultihasenCancleMonomialParams(xml, path));
				}
			bool registered = false;
		}

		const std::string name("TWO_FLAVOR_SEOPREC_CONSTDET_MULTIHASEN_CANCLE_FERM_MONOMIAL");
		bool registerAll()
		{
			bool success = true;
			if(!registered)
			{
				success &= WilsonTypeFermActs4DEnv::registerAll();
				success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
				registered = true;
			}
			return success;
		}
	}
	// Constructor
	SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial::
		SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(const
				TwoFlavorMultihasenCancleMonomialParams& param)
		{
			START_CODE();

			inv_param = param.inv_param;
			std::istringstream is(param.fermact.xml);
			XMLReader fermact_reader(is);
			QDPIO::cout<<__func__<<": construct "
				<<param.fermact.id<<std::endl;

			WilsonTypeFermAct<T,P,Q>* tmp_act = 
				TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id,fermact_reader, 
						param.fermact.path);
			ShiftSymEvenOddPrecCloverFermAct* downcast=dynamic_cast<ShiftSymEvenOddPrecCloverFermAct*>(tmp_act);

			if(downcast == 0x0)
			{
				QDPIO::cerr<<__func__<<": unable to downcast FermAct to ShiftSymEvenOddPrecCloverFermAct"<<std::endl;
				QDP_abort(1);
			}
			fermact = downcast;

			mu = param.mu;
			// Set shifted mass parameter for ferm action
			fermact->setMu(mu);

			// Get Chronological predictor
			AbsChronologicalPredictor4D<LatticeFermion>* tmp = 0x0;
			if(param.predictor.xml == ""){
				tmp = new ZeroGuess4DChronoPredictor();
			}else{
				try{
					std::istringstream chrono_is(param.predictor.xml);
					XMLReader chrono_xml(chrono_is);
					tmp = The4DChronologicalPredictorFactory::Instance().createObject(param.predictor.id,
							chrono_xml, param.predictor.path);
				}catch(const std::string& e){
					QDPIO::cerr<<"Caught Exception Reading XML: "<<e<<std::endl;
					QDP_abort(1);
				}
			}
			if(tmp==0x0){
				QDPIO::cerr<<"Failed to create ZeroGuess4DChronoPredictor"<<std::endl;
				QDP_abort(1);
			}
			chrono_predictor = tmp;
			END_CODE();
		}

	Double SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial::
		S_odd_odd(const AbsFieldState<multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >& s){
			START_CODE();

			XMLWriter& xml_out = TheXMLLogWriter::Instance();
			push(xml_out, "S_odd_odd");

			// Fermion action
			ShiftSymEvenOddPrecCloverFermAct& FA = getFermAct();

			Handle<FermState<T,P,Q> > state = FA.createState(s.getQ());

			// Linear Operator
			Handle<SymEvenOddPrecLogDetLinearOperator<T,P,Q> > lin(FA.linOp(state));

			// Get system solver
			Handle<MdagMSystemSolver<T> > invMdagM(FA.invMdagM(state, getInvParams()));

			// Get the X fields
			T X;
			X[lin->subset()] = zero;

			// No predictor used here.
			QDPIO::cout<<"TwoFlavMultihasenCancleMonomial: resetting Predictor before energy calc solve" << std::endl;
			(getMDSolutionPredictor()).reset();

			SystemSolverResults_t res = (*invMdagM)(X, getPhi());
			QDPIO::cout<<"2Flav Multihasen Cancle, n_count = "<< res.n_count<<std::endl;

			LatticeDouble site_action = zero;
			site_action[lin->subset()] = Double(-12);
			site_action[lin->subset()] += localInnerProductReal(getPhi(), X);

			Double action = sum(site_action, lin->subset());

			write(xml_out, "n_count", res.n_count);
			write(xml_out, "S_oo", action);
			pop(xml_out);

			END_CODE();

			return action;
		}

	void SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial::
		dsdq(P& F, const AbsFieldState<P,Q>& s)
		{
			START_CODE();

			XMLWriter& xml_out = TheXMLLogWriter::Instance();
			push(xml_out, "TwoFlavorConstDetMultihasenCancleMonomial");

			ShiftSymEvenOddPrecCloverFermAct& FA = getFermAct();
			Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));

			Handle<MdagMSystemSolver<T> > invMdagM(FA.invMdagM(state, getInvParams()));
			Handle<SymEvenOddPrecLogDetLinearOperator<T,P,Q> > M(FA.linOp(state));

			T X;
			SystemSolverResults_t res = (*invMdagM)(X, getPhi(), getMDSolutionPredictor());
			QDPIO::cout<<"2FlavConstDetMultihasenCancle:: invert, n_count = "<<res.n_count<<std::endl;
		
			T Y;
			(*M)(Y, X, PLUS);
			M->deriv(F, X, Y, MINUS);

			P F_tmp;
			M->deriv(F_tmp, Y, X, PLUS);
			F += F_tmp;

			for(int i=0; i<F.size(); ++i)
				F[i] *= Real(-1);
			
			state->deriv(F);
			
			write(xml_out, "n_count", res.n_count);
			monitorForces(xml_out, "Forces", F);
			pop(xml_out);

			END_CODE();
		}

	void SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial::
		refreshInternalFields(const AbsFieldState<P,Q>& field_state)
		{
			START_CODE();
			ShiftSymEvenOddPrecCloverFermAct& FA = getFermAct();
			Handle<FermState<T,P,Q> > state(FA.createState(field_state.getQ()));
			Handle<SymEvenOddPrecLogDetLinearOperator<T,P,Q> > M(FA.linOp(state));

			T eta = zero;
			gaussian(eta, M->subset());
			FA.getFermBC().modifyF(eta);

			eta *= sqrt(0.5);
			(*M)(getPhi(), eta, MINUS);

			QDPIO::cout<<__func__<<": resettng predictor after field refresh"<<std::endl;
			getMDSolutionPredictor().reset();

			END_CODE();
		}

	void SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial::
		setInternalFields(const Monomial<P,Q>& m)
		{
			START_CODE();
			try{
				const SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial& fm = 
					dynamic_cast<const SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial&>(m);
				getPhi() = fm.getPhi();
			}catch(std::bad_cast){
				QDPIO::cerr<<"Failed to cast input Monomial to SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial"<<
					std::endl;
				QDP_abort(1);
			}

			END_CODE();
		}

}


