/* Multi-Hasenbusch cancle monomial term for both 
 * symm and asymm even-odd precondition
 */
#ifndef __TWO_FLAVOR_MULTIHASEN_CANCLE_MONOMIAL_W_H__
#define __TWO_FLAVOR_MULTIHASEN_CANCLE_MONOMIAL_W_H__

#include "update/molecdyn/field_state.h"
#include "actions/ferm/linop/shifted_linop_w.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "chromabase.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "seoprec_constdet_wilstype_fermact_w.h"
#include "seoprec_logdet_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "eoprec_logdet_wilstype_fermact_w.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/monomial/two_flavor_multihasen_cancle_monomial_params_w.h"


namespace Chroma
{
	// Symmetric preconditioned 
	namespace SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomialEnv
	{
		bool registerAll();
	}

	// Asymmetric preconditioned
	namespace EvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomialEnv
	{
		bool registerAll();
	}

	template<typename T, typename P, typename Q,
		template<typename, typename, typename> class FAType,
		template<typename, typename, typename> class LOType>
			class PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial: 
				public ExactWilsonTypeFermMonomial<P,Q,T>
	{
		// same as old monomial body but replace the specific FermionAction
		// typename with "FAType" and every specific LinearOperator type with
		// LOType
		public:
			PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(const
					TwoFlavorMultihasenCancleMonomialParams& param_);
			PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(const
					PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial& m):
				phi(m.phi), fermact(m.fermact), inv_param(m.inv_param),
				chrono_predictor(m.chrono_predictor){}

			virtual ~PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(){}

			virtual Double S(const AbsFieldState<P,Q>& s)
			{
				START_CODE();

				XMLWriter& xml_out = TheXMLLogWriter::Instance();
				push(xml_out, "S");

				const FAType<T,P,Q>& FA = getFermAct();
				Handle<FermState<T,P,Q> > state = FA.createState(s.getQ());

				// Get the X fields
				T X;
				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));

				// Shifted lineart operator with mass parameter mu
				Handle<LinearOperator<T> >
					M(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu));
				X[M->subset()] = zero;
				// Energy calc doesnt use Chrono Predictor
				QDPIO::cout << "TwoFlavWilson4DCancleMonomial: resetting Predictor before energy calc solve" << std::endl;
				(getMDSolutionPredictor()).reset();

				// Get system solver
				const GroupXML_t& invParam = getInvParams();
				std::istringstream xml(invParam.xml);
				XMLReader paramtop(xml);
				Handle<MdagMSystemSolver<T> >
					invMdagM(TheMdagMFermSystemSolverFactory::Instance().createObject(
								invParam.id, paramtop, invParam.path, state, M));
				// Solve MdagM X = eta
				SystemSolverResults_t res = (*invMdagM)(X, getPhi());
				QDPIO::cout<<"2FlavCancle::invert, n_count = "<<res.n_count<<std::endl;

				// Action
				Double action = innerProductReal(getPhi(), X, M->subset());
				write(xml_out, "n_count", res.n_count);
				write(xml_out, "S", action);
				pop(xml_out);

				END_CODE();
				return action;
			}

			virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
			{
				START_CODE();

				XMLWriter& xml_out = TheXMLLogWriter::Instance();
				push(xml_out, "TwoFlavorWilsonTypeMultihasenCancleMonomial");

				const FAType<T,P,Q>& FA = getFermAct();
				Handle<FermState<T,P,Q> > state = FA.createState(s.getQ());

				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));

				// Shifted lineart operator with mass parameter mu
				Handle<LinearOperator<T> >
					M(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu));

				// Get system solver
				const GroupXML_t& invParam = getInvParams();
				std::istringstream xml(invParam.xml);
				XMLReader paramtop(xml);
				Handle<MdagMSystemSolver<T> >
					invMdagM(TheMdagMFermSystemSolverFactory::Instance().createObject(
								invParam.id, paramtop, invParam.path, state, M));
				T X;
				// Solve MdagM X = eta
				SystemSolverResults_t res = (*invMdagM)(X, getPhi(),getMDSolutionPredictor());
				QDPIO::cout << "2FlavCancle::invert,  n_count = " << res.n_count << std::endl;

				P F_tmp;
				T Y;
				(*M)(Y, X, PLUS);

				// cast linearop to difflinearop
				const DiffLinearOperator<T,P,Q>& diffM = 
					dynamic_cast<const DiffLinearOperator<T,P,Q>&>(*M);
				diffM.deriv(F, X, Y, MINUS);
				diffM.deriv(F_tmp, Y, X, PLUS);
				F += F_tmp;

				for(int i=0; i<Nd; ++i){
					F[i] *= Real(-1);
				}

				state->deriv(F);

				write(xml_out, "n_count", res.n_count);
				monitorForces(xml_out, "Forces", F);

				pop(xml_out);
				END_CODE();
			}

			virtual void refreshInternalFields(const AbsFieldState<P,Q>& s)
			{
				START_CODE();

				const FAType<T,P,Q>& FA = getFermAct();
				Handle<FermState<T,P,Q> > state = FA.createState(s.getQ());

				Handle<LOType<T,P,Q> > base_op(FA.linOp(state));

				// Shifted lineart operator with mass parameter mu
				Handle<LinearOperator<T> >
					M(new TwistedShiftedLinOp<T,P,Q,LOType>(*base_op, mu));

				T eta = zero;
				gaussian(eta, M->subset());
				FA.getFermBC().modifyF(eta);

				eta *= sqrt(0.5);
				(*M)(getPhi(), eta, MINUS);

				QDPIO::cout << "TwoFlavWilson4DCancleMonomial: resetting Predictor after field refresh" << std::endl;
				getMDSolutionPredictor().reset();

				END_CODE();
			}

			virtual void setInternalFields(const Monomial<P,Q>& m)
			{
				START_CODE();
				try{
					const PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial<T,P,Q,FAType,LOType>& fm = 
						dynamic_cast<const PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial<T,P,Q,FAType,LOType>& >(m);
					getPhi() = fm.getPhi();
				}
				catch(std::bad_cast){
					QDPIO::cerr<<"Failed to cast input Monomial to PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial "<<std::endl;
					QDP_abort(1);
				}
				END_CODE();
			}
			virtual void resetPredictors(void){
				getMDSolutionPredictor().reset();
			}

		protected:
			virtual const T& getPhi(void) const {
				return phi;
			}
			virtual T& getPhi(void){
				return phi;
			}
			const FAType<T,P,Q>& getFermAct(void) const{
				return *fermact;
			}
			const GroupXML_t& getInvParams(void) const{
				return inv_param;
			}
			AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void){
				return *chrono_predictor;
			}

		private:
			PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial();
			void operator=(const PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial&);
			// Pseudofermion field phi
			T phi;
			Handle<const FAType<T,P,Q> > fermact;
			// Shifted mass parameter
			Real mu;
			GroupXML_t inv_param;
			Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
	};

	template<typename T, typename P, typename Q,
		template<typename, typename, typename> class FAType,
		template<typename, typename, typename> class LOType> 
			PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial<T,P,Q,FAType,LOType>::
			PrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(const
					TwoFlavorMultihasenCancleMonomialParams& param)
			{
				START_CODE();
				inv_param = param.inv_param;
				std::istringstream is(param.fermact.xml);
				XMLReader fermact_reader(is);
				QDPIO::cout<<__func__<<": construct "
					<<param.fermact.id<<std::endl;

				WilsonTypeFermAct<T,P,Q>* tmp_act = 
					TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id,
							fermact_reader,
							param.fermact.path);
				FAType<T,P,Q>* downcast=dynamic_cast<FAType<T,P,Q>* >(tmp_act);

				if(downcast == 0x0){
					QDPIO::cerr<<__func__<<": unable to downcast FermAct"<<std::endl;
					QDP_abort(1);
				}
				fermact = downcast;

				mu = param.mu;

				AbsChronologicalPredictor4D<LatticeFermion>* tmp = 0x0;
				if(param.predictor.xml == ""){
					tmp = new ZeroGuess4DChronoPredictor();
				}else{
					try{
						std::istringstream chrono_is(param.predictor.xml);
						XMLReader chrono_xml(chrono_is);
						tmp = The4DChronologicalPredictorFactory::Instance().createObject(param.predictor.id,
								chrono_xml,
								param.predictor.path);
					}catch(const std::string& e){
						QDPIO::cerr<<"Caught Excepthion Reading XML: "<<e<<std::endl;
						QDP_abort(1);
					}
				}
				if(tmp == 0x0){
					QDPIO::cerr<<"Failed to create ZeroGuess4DChronoPredictor"<<std::endl;
					QDP_abort(1);
				}
				chrono_predictor = tmp;
				END_CODE();
			}

}//end namespace Chroma

#endif
