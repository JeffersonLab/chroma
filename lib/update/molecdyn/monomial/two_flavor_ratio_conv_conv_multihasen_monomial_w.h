// -*- C++ -*-
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_W_H__
#define __TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_W_H__

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_multihasen_monomial_params_w.h"
#include "actions/ferm/linop/shifted_linop_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "seoprec_constdet_wilstype_fermact_w.h"
#include "seoprec_logdet_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "eoprec_logdet_wilstype_fermact_w.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"

namespace Chroma
{

	// Symmetric preconditioned
	namespace SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialEnv
	{
		bool registerAll();
	}

	// Asymmetric preconditioned
	namespace EvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialEnv
	{
		bool registerAll();
	}

	template<typename T, typename P, typename Q,
		template<typename, typename, typename> class FAType,
		template<typename, typename, typename> class LOType>
			class PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial: 
				public ExactWilsonTypeFermMonomial<P,Q,T>
	{
		public:
			PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(const 
					TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param_);
			virtual ~PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(){}

			virtual Double S(const AbsFieldState<P, Q>& s)
			{
				START_CODE();
				XMLWriter& xml_out = TheXMLLogWriter::Instance();
				push(xml_out, "PrecConstDetTwoFlavorRatioConvConvMultihasenFermMonomial");

				// Fermion action
				const FAType<T,P,Q>& FA = getFermAct();

				Double S=0;
				// Get the X fields
				T X;
				Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));

				for(int i=0; i<numHasenTerms; ++i){
					Handle<LinearOperator<T> > 
						M(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu[i]));
					Handle<LinearOperator<T> > 
						M_prec(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu[i+1]));

					// Action calc doesnt use chrono predictor use zero guess
					X[M->subset()] = zero;
					// Reset chrono predictor
					QDPIO::cout << 
						"TwoFlavorRatioConvConvMultihasenWilson4DMonomial: resetting Predictor before energy calc solve" << std::endl;
					(getMDSolutionPredictor()).reset();
					// Get system solver
					const GroupXML_t& invParam = getInvParams();
					std::istringstream xml(invParam.xml);
					XMLReader paramtop(xml);
					Handle<MdagMSystemSolver<T> > 
						invMdagM(TheMdagMFermSystemSolverFactory::Instance().createObject(
									invParam.id, paramtop, invParam.path, state, M));
					// M_dag_prec phi = M^{dag}_prec \phi - the RHS
					T M_dag_prec_phi;
					(*M_prec)(M_dag_prec_phi, getPhi(i), MINUS);

					// Solve MdagM X = eta
					SystemSolverResults_t res = (*invMdagM)(X, M_dag_prec_phi);

					T phi_tmp = zero;
					(*M_prec)(phi_tmp, X, PLUS);
					Double action = innerProductReal(getPhi(i), phi_tmp, M->subset());
					// Write out inversion number and action for every hasenbusch term 
					std::string n_count = "n_count_hasenterm_" + std::to_string(i+1);
					std::string s_action = "S_hasenterm_" + std::to_string(i+1);

					write(xml_out, n_count, res.n_count);
					write(xml_out, s_action, action);
					S += action;
				}
				write(xml_out, "S", S);
				pop(xml_out);

				END_CODE();
				return S;
			}

			// Total force of multi-hasen terms
			virtual void dsdq(P& F, const AbsFieldState<P, Q>& s)
			{
				START_CODE();
				XMLWriter& xml_out = TheXMLLogWriter::Instance();
				push(xml_out, "TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial");

				P F_t;
				P F_tmp;
				F_t.resize(Nd);
				F_tmp.resize(Nd);
				F.resize(Nd);
				for(int i=0; i<Nd; ++i)
				{
					F_t[i] = zero;
					F_tmp[i] = zero;
					F[i] = zero;
				}
				// Fermion action
				const FAType<T,P,Q>& FA = getFermAct();
				// X = (M^\dag M)^(-1) M_prec^\dag \phi
				// Y = M X
				T X = zero;
				T Y = zero;
				// M_dag_prec_phi = M^\dag_prec \phi
				T M_dag_prec_phi;

				Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));

				for(int i=0; i<numHasenTerms; ++i){
					Handle<LinearOperator<T> > 
						M(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu[i]));
					Handle<LinearOperator<T> > 
						M_prec(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu[i+1]));
					// Get system solver
					const GroupXML_t& invParam = getInvParams();
					std::istringstream xml(invParam.xml);
					XMLReader paramtop(xml);
					Handle<MdagMSystemSolver<T> > 
						invMdagM(TheMdagMFermSystemSolverFactory::Instance().createObject(
									invParam.id,paramtop, invParam.path, state, M));

					(*M_prec)(M_dag_prec_phi, getPhi(i), MINUS);
					SystemSolverResults_t res = (*invMdagM)(X, M_dag_prec_phi, getMDSolutionPredictor());
					(*M)(Y, X, PLUS);

					// cast linearop to difflinearop
					const DiffLinearOperator<T,P,Q>& diffM = dynamic_cast<const DiffLinearOperator<T,P,Q>&>(*M);
					const DiffLinearOperator<T,P,Q>& diffM_prec = dynamic_cast<const DiffLinearOperator<T,P,Q>&>(*M_prec);

					// deriv part 1: \phi^\dag \dot(M_prec) X
					diffM_prec.deriv(F_t, getPhi(i), X, PLUS);

					// deriv part 2: - X^\dag \dot(M^\dag) Y
					diffM.deriv(F_tmp, X, Y, MINUS);
					F_t -= F_tmp;

					// deriv part 3: - Y^\dag \dot(M) X
					diffM.deriv(F_tmp, Y, X, PLUS);
					F_t -= F_tmp;

					// deriv part 4: X^\dag \dot(M_prec)^\dag \phi
					diffM_prec.deriv(F_tmp, X, getPhi(i), MINUS);
					F_t += F_tmp;

					// total force from all Hasenbusch terms
					F += F_t;

					// F now holds derivative with respect to possibly fat links
					// now derive it with respect to the thin links if needs be
					state->deriv(F);

					// Write out inversion number and action for every hasenbusch term 
					std::string n_count = "n_count_hasenterm_" + std::to_string(i+1);
					std::string force = "Forces_hasenterm_" + std::to_string(i+1);

					write(xml_out, n_count, res.n_count);
					monitorForces(xml_out, force, F_t);
				}
				// Total force from all Hasenbusch terms
				monitorForces(xml_out, "Forces", F);
				pop(xml_out);
				END_CODE();
			}

			virtual void refreshInternalFields(const AbsFieldState<P, Q>& s)
			{
				START_CODE();
				const FAType<T,P,Q>& FA = getFermAct();

				Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));
				for(int i=0; i<numHasenTerms; ++i){
					Handle<LinearOperator<T> > 
						M(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu[i]));
					Handle<LinearOperator<T> > 
						M_prec(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu[i+1]));
					T eta = zero;
					gaussian(eta, M->subset());

					// Account for fermion BC by modifying the proposed field
					FA.getFermBC().modifyF(eta);

					eta *= sqrt(0.5);
					// Now we have to get the refreshment right:
					// We have  \eta^{\dag} \eta = \phi M_prec (M^dag M )^-1 M^dag_prec \phi
					//  so that \eta = (M^\dag)^{-1} M^{dag}_prec \phi
					//  So need to solve M^{dag}_prec \phi = M^{\dag} \eta
					// Which we can solve as 
					//      M^{dag}_prec M_prec ( M_prec^{-1} ) \phi = M^{\dag} \eta
					//
					// Zero out everything - to make sure there is no junk
					// in uninitialised places
					T eta_tmp = zero;
					T phi_tmp = zero;
					getPhi(i) = zero;

					(*M)(eta_tmp, eta, MINUS); // M^\dag \eta

					// Get system solver
					const GroupXML_t& invParam = getInvParams();
					std::istringstream xml(invParam.xml);
					XMLReader paramtop(xml);
					Handle<MdagMSystemSolver<T> > 
						invMdagM(TheMdagMFermSystemSolverFactory::Instance().createObject(
									invParam.id,paramtop, invParam.path, state, M_prec));

					// Solve MdagM_prec X = eta
					SystemSolverResults_t res = (*invMdagM)(phi_tmp, eta_tmp);
					(*M_prec)(getPhi(i), phi_tmp, PLUS); // (Now get phi = M_prec (M_prec^{-1}\phi)
					QDPIO::cout<<"TwoFlavRatioConvConvMultihasenWilson4DMonomial: resetting Predictor at end of field refresh"<<std::endl;
					getMDSolutionPredictor().reset();
					XMLWriter& xml_out = TheXMLLogWriter::Instance();
					std::string field_refrs = "FieldRefreshment_"+std::to_string(i);
					std::string n_count = "n_count_"+std::to_string(i);
					push(xml_out, field_refrs);
					write(xml_out, n_count, res.n_count);
					pop(xml_out);
				}
				END_CODE();
			}

			virtual void setInternalFields(const Monomial<P, Q>& m)
			{
				START_CODE();
				try{
					const PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial<T,P,Q,FAType,LOType>& fm =
						dynamic_cast<const PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial<T,P,Q,FAType,LOType>& >(m);
					for(int i=0; i<numHasenTerms; ++i)
						getPhi(i) = fm.getPhi(i);
				}catch(std::bad_cast){
					QDPIO::cerr<<"Failed to cast input Monomial to PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial"
						<<std::endl;
					QDP_abort(1);
				}
				END_CODE();
			}

			void resetPredictors(){
				getMDSolutionPredictor().reset();
			}

			int getNumHasenTerms() const {
				return numHasenTerms;
			}

		protected:
			// override base class getPhi()
			virtual const T& getPhi(void) const {};
			virtual T& getPhi(void) {};
			// pesudo-fermion field for each Hasenbusch term
			virtual T& getPhi(int i) {
				return phi[i];
			}
			virtual const T& getPhi(int i) const {
				return phi[i];
			}

			const FAType<T,P,Q>& getFermAct()const{
				return *fermact;
			}

			//! Get parameters for the inverter
			const GroupXML_t& getInvParams() const { 
				return invParam;
			}

			AbsChronologicalPredictor4D<T>& getMDSolutionPredictor() { 
				return *chrono_predictor;
			}
		private:
			// Hide empty constructor and =
			PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial();
			void operator=(const PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial&);

			// Pseudofermion field phi for multi-hasenbusch term
			multi1d<T> phi;

			Handle<const FAType<T,P,Q> > fermact;

			// Shifted mass mu
			multi1d<Real> mu;

			// Number of Hasenbusch terms
			int numHasenTerms;

			// The parameters for the inversion
			GroupXML_t invParam;

			Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
	};


	template<typename T, typename P, typename Q,
		template<typename, typename, typename> class FAType,
		template<typename, typename, typename> class LOType> 
			PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial<T,P,Q,FAType,LOType>::
			PrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(const
					TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param)
			{
				START_CODE();
				QDPIO::cout << "Constructor: " << __func__ << std::endl;

				if( param.fermactInv.invParam.id == "NULL" ) {
					QDPIO::cerr << "Fermact inverter params are NULL" << std::endl;
					QDP_abort(1);
				}

				invParam = param.fermactInv.invParam;

				// Fermion action (with original seoprec action)
				{
					std::istringstream is(param.fermactInv.fermact.xml);
					XMLReader fermact_reader(is);
					QDPIO::cout << "Construct fermion action= " << param.fermactInv.fermact.id << std::endl;

					WilsonTypeFermAct<T,P,Q>* tmp_act = 
						TheWilsonTypeFermActFactory::Instance().createObject(param.fermactInv.fermact.id, 
								fermact_reader, 
								param.fermactInv.fermact.path);

					FAType<T,P,Q>* downcast=dynamic_cast<FAType<T,P,Q>* >(tmp_act);

					// Check success of the downcast 
					if( downcast == 0x0 ) 
					{
						QDPIO::cerr << __func__ << ": unable to downcast FermAct" << std::endl;
						QDP_abort(1);
					}
					fermact = downcast;    
				}

				// Number of Hasenbusch term
				numHasenTerms = param.numHasenTerms;

				// Shifted mass
				mu.resize(numHasenTerms+1);
				for(int i=0; i< numHasenTerms+1; ++i){
					mu[i] = param.mu[i];
				}

				// Pseudofermion field phi
				phi.resize(numHasenTerms);
				for(int i=0; i<numHasenTerms; ++i)
					phi[i] = zero;

				// Get Chronological predictor
				{
					AbsChronologicalPredictor4D<LatticeFermion>* tmp = 0x0;
					if( param.predictor.xml == "" ) {
						// No predictor specified use zero guess
						tmp = new ZeroGuess4DChronoPredictor();
					}else {
						try { 
							std::istringstream chrono_is(param.predictor.xml);
							XMLReader chrono_xml(chrono_is);
							tmp = The4DChronologicalPredictorFactory::Instance().createObject(param.predictor.id, 
									chrono_xml, 
									param.predictor.path);
						}catch(const std::string& e ) { 
							QDPIO::cerr << "Caught Exception Reading XML: " << e << std::endl;
							QDP_abort(1);
						}
					}
					if( tmp == 0x0 ) { 
						QDPIO::cerr << "Failed to create the 4D ChronoPredictor" << std::endl;
						QDP_abort(1);
					}
					chrono_predictor = tmp;
				}

				QDPIO::cout << "Finished constructing: " << __func__ << std::endl;
				END_CODE();
			}

} //end namespace chroma
#endif
