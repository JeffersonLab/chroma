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
#include "actions/ferm/invert/syssolver_mrhs_twisted_params.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_proxy.h"
#include "actions/ferm/invert/syssolver_linop_mrhs_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_mrhs_factory.h"
#include "actions/ferm/linop/multi_twist_linop_w.h"

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

				const int N = numHasenTerms;
				Double S=0;
				// Get the X fields
				multi1d<T> X(N);
				Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));
				
				// Construct Linop with multi-twist
				// Denominator of multi-Hasenbusch monomial
				Handle<MultiTwistLinOp<T,P,Q,LOType> >
					M_den(new MultiTwistLinOp<T,P,Q,LOType>(base_op, mu_den));
				// Numerator of multi-Hasenbusch monomial, not needed here,
				// use fermion action directly.
				//Handle<MultiTwistLinOp<T,P,Q,LOType> >
				//	M_num(new MultiTwistLinOp<T,P,Q,LOType>(base_op, mu_num));
				
				// Action calc doesnt use chrono predictor use zero guess
				for(int i=0; i<N; ++i){
					X[i][M_den->subset()] = zero;
				}
				// Reset chrono predictor
				QDPIO::cout << 
					"TwoFlavorRatioConvConvMultihasenWilson4DMonomial: resetting Predictor before energy calc solve" << std::endl;
				(getMDSolutionPredictor()).reset();
			
				// A bit tedious work to set up the Twists for numerator
				// MdagM solver
				SysSolverMRHSTwistedParams paramnum = invParam;
				paramnum.Twists = mu_num;
				// Set up the MRHS solver for MdagM, where M is the
				// twisted linop on numerator.
				MdagMMRHSSysSolverTwistedProxy<T,P,Q,LOType> invMdagM(paramnum, FA, state);
				
				// M_dag_den_phi = M^{dag}_den \phi - the RHS
				multi1d<T> M_dag_den_phi(N);
				(*M_den)(M_dag_den_phi, phi, MINUS);

				// Solve MdagM X = eta
				SystemSolverResultsMRHS_t res = invMdagM(X, M_dag_den_phi);
				multi1d<T> phi_tmp(N);
				phi_tmp = zero;
				(*M_den)(phi_tmp, X, PLUS);

				for(int i=0; i<N; ++i){

					Double action = innerProductReal(getPhi(i), phi_tmp[i], M_den->subset());
					// Write out inversion number and action for every hasenbusch term 
					std::string n_count = "n_count_hasenterm_" + std::to_string(i+1);
					std::string s_action = "S_hasenterm_" + std::to_string(i+1);

					write(xml_out, n_count, res.n_count[i]);
					write(xml_out, s_action, action);
					S += action;
				}
				write(xml_out, "S", S);
				pop(xml_out);

				END_CODE();
				return S;
			}

			// Total force of multi-hasen terms
			virtual void dsdq(P& F_t, const AbsFieldState<P, Q>& s)
			{
				START_CODE();
				XMLWriter& xml_out = TheXMLLogWriter::Instance();
				push(xml_out, "TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial");

				F_t.resize(Nd);
				for(int i=0; i<Nd; ++i)
				{
					F_t[i] = zero;
				}
				const int N = numHasenTerms;
				multi1d<P> F_tmp1(N), F_tmp2(N);
                for(int i=0; i<N; ++i){
                    F_tmp1[i].resize(Nd);
                    F_tmp2[i].resize(Nd);
                    for(int j=0; j<Nd; ++j){
                        F_tmp1[i][j] = zero;
                        F_tmp2[i][j] = zero;
                    }
                }
                // Fermion action
				const FAType<T,P,Q>& FA = getFermAct();
				// X = (M^\dag M)^(-1) M_den^\dag \phi
				// Y = M X
				multi1d<T> X(N), Y(N);
				for(int i=0; i<N; ++i){
					X[i] = zero;
					Y[i] = zero;
				}
                // M_dag_den_phi = M^\dag_den \phi
				multi1d<T> M_dag_den_phi(N);

				Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));

				// Construct Linop with multi-twist
				// Denominator of multi-Hasenbusch monomial
				Handle<MultiTwistLinOp<T,P,Q,LOType> >
					M_den(new MultiTwistLinOp<T,P,Q,LOType>(base_op, mu_den));
				// Numerator of multi-Hasenbusch monomial
				Handle<MultiTwistLinOp<T,P,Q,LOType> >
					M_num(new MultiTwistLinOp<T,P,Q,LOType>(base_op, mu_num));
				
				// A bit tedious work to set up the Twists for numerator
				// MdagM solver
				SysSolverMRHSTwistedParams paramnum = invParam;
				paramnum.Twists = mu_num;
				// Set up the MRHS solver for MdagM, where M is the
				// twisted linop on numerator.
				MdagMMRHSSysSolverTwistedProxy<T,P,Q,LOType> invMdagM(paramnum, FA, state);
			
				(*M_den)(M_dag_den_phi, phi, MINUS);
				SystemSolverResultsMRHS_t res = invMdagM(X, M_dag_den_phi);
				(*M_num)(Y, X, PLUS);
				
				// deriv part 1: \phi^\dag \dot(M_den) X
				M_den->deriv(F_tmp1, phi, X, PLUS);

				// deriv part 2: - X^\dag \dot(M_num^\dag) Y
				M_num->deriv(F_tmp2, X, Y, MINUS);
				F_tmp1 -= F_tmp2;

				// deriv part 3: - Y^\dag \dot(M_num) X
				M_num->deriv(F_tmp2, Y, X, PLUS);
				F_tmp1 -= F_tmp2;

				// deriv part 4: X^\dag \dot(M_prec)^\dag \phi
				M_den->deriv(F_tmp2, X, phi, MINUS);
				F_tmp1 += F_tmp2;
    
                // Total force
                for(int i=0; i<N; ++i){
                    F_t += F_tmp1[i];
                }
				// F now holds derivative with respect to possibly fat links
				// now derive it with respect to the thin links if needs be
				state->deriv(F_t);

				for(int i=0; i<N; ++i){
					// Write out inversion number for every hasenbusch term 
					std::string n_count = "n_count_hasenterm_" + std::to_string(i+1);
					write(xml_out, n_count, res.n_count[i]);
                    monitorForces(xml_out, "Forces_hasenterm_" + std::to_string(i+1), F_tmp1[i]);
				}
				// Total force from all Hasenbusch terms
				monitorForces(xml_out, "ForcesTotal", F_t);
				pop(xml_out);
				END_CODE();
			}

			virtual void refreshInternalFields(const AbsFieldState<P, Q>& s)
			{
				START_CODE();
				const FAType<T,P,Q>& FA = getFermAct();
				const int N = numHasenTerms;
				
				Handle<FermState<T,P,Q> > state(FA.createState(s.getQ()));
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));
				// Construct Linop with multi-twist
				// Denominator of multi-Hasenbusch monomial
				Handle<MultiTwistLinOp<T,P,Q,LOType> >
					M_den(new MultiTwistLinOp<T,P,Q,LOType>(base_op, mu_den));
				// Numerator of multi-Hasenbusch monomial
				Handle<MultiTwistLinOp<T,P,Q,LOType> >
					M_num(new MultiTwistLinOp<T,P,Q,LOType>(base_op, mu_num));
				
				multi1d<T> eta(N);
				for(int i=0; i<N; ++i){
					eta[i] = zero;
					gaussian(eta[i], M_num->subset());
					FA.getFermBC().modifyF(eta[i]);
					eta[i] *= sqrt(0.5);
				}
				// Now we have to get the refreshment right:
				// We have  \eta^{\dag} \eta = \phi^{\dag} M_den (M_num^dag M_num )^-1 M^dag_den \phi
				//  so that \eta = (M_num^\dag)^{-1} M^{dag}_den \phi
				//  So need to solve M^{dag}_den \phi = M_num^{\dag} \eta
				// Which we can solve as 
				//      M^{dag}_den M_den ( M_den^{-1} ) \phi = M_num^{\dag} \eta
				//
				// Zero out everything - to make sure there is no junk
				// in uninitialised places
				multi1d<T> eta_tmp(N), phi_tmp(N);
				for(int i=0; i<N; ++i){
					eta_tmp[i] = zero;
					phi_tmp[i] = zero;
					getPhi(i) = zero;
				}
				(*M_num)(eta_tmp, eta, MINUS); // M_num^\dag \eta
			
				// Get system solver
				//SysSolverMRHSTwistedParams paramnum = invParam;
				//paramnum.Twists = mu_num;
				// Set up the MRHS solver for MdagM, here for denomerator.
				MdagMMRHSSysSolverTwistedProxy<T,P,Q,LOType> invMdagM(invParam, FA, state);
				// Solve MdagM_den X = eta
				SystemSolverResultsMRHS_t res = invMdagM(phi_tmp, eta_tmp);
				(*M_den)(phi, phi_tmp, PLUS); // (Now get phi = M_den (M_den^{-1}\phi)
				QDPIO::cout<<"TwoFlavRatioConvConvMultihasenWilson4DMonomial: resetting Predictor at end of field refresh"<<std::endl;
				getMDSolutionPredictor().reset();
				
				XMLWriter& xml_out = TheXMLLogWriter::Instance();
				for(int i=0; i<N; ++i){
					std::string field_refrs = "FieldRefreshment_"+std::to_string(i);
					std::string n_count = "n_count_"+std::to_string(i);
					push(xml_out, field_refrs);
					write(xml_out, n_count, res.n_count[i]);
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

			// Shifted mass from denominator of multi-Hasenbusch term
			multi1d<Real> mu_den;
			// Shifted mass from numerator of multi-Hasenbusch term
			multi1d<Real> mu_num;
			// No. of Hasenbusch terms
			int numHasenTerms;
			// The parameters for the inversion
			SysSolverMRHSTwistedParams invParam;

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

				if( param.fermact.id == "NULL" ) {
					QDPIO::cerr << "Fermact inverter params are NULL" << std::endl;
					QDP_abort(1);
				}

				invParam = param.inverter;

				// Fermion action (with original seoprec action)
				{
					std::istringstream is(param.fermact.xml);
					XMLReader fermact_reader(is);
					QDPIO::cout << "Construct fermion action= " << param.fermact.id << std::endl;

					WilsonTypeFermAct<T,P,Q>* tmp_act = 
						TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id, 
								fermact_reader, 
								param.fermact.path);

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
				numHasenTerms = param.inverter.BlockSize;
				// Shifted mass of denominator
				mu_den.resize(numHasenTerms);
				for(int i=0; i< numHasenTerms; ++i){
					mu_den[i] = param.inverter.Twists[i];
				}
				// Shifted mass of numerator
				mu_num.resize(numHasenTerms);
				mu_num[0] = param.base_twist;
				for(int i=1; i< numHasenTerms; ++i){
					mu_num[i] = param.inverter.Twists[i-1];
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
