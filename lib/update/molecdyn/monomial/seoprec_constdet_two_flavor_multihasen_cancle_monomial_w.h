// -*- C++ -*-
/*! @file
 * @brief Two-flavor collection of symmetric even-odd preconditioned 4D ferm monomials
 * cancle term for multi-Hasenbusch ratio 
 */

#ifndef __SEOPREC_CONSTDET_TWO_FLAVOR_MULTIHASEN_CANCLE_H_
#define __SEOPREC_CONSTDET_TWO_FLAVOR_MULTIHASEN_CANCLE_H_

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "actions/ferm/fermacts/shifted_seoprec_clover_fermact_w.h"
#include "update/molecdyn/monomial/two_flavor_multihasen_cancle_monomial_params_w.h"

namespace Chroma
{
	namespace SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomialEnv
	{
		bool registerAll();
	}

	// Inherit and utilize SymEvenOddPrecConstDet...Monomial class
	// only difference is construction of fermact with shifted mass mu
	class SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial:
		public TwoFlavorExactSymEvenOddPrecConstDetWilsonTypeFermMonomial<
		multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix>,
		LatticeFermion>
	{
		public:
			typedef LatticeFermion				 T;
			typedef multi1d<LatticeColorMatrix>  P;
			typedef multi1d<LatticeColorMatrix>  Q;
			
			// Construct out of a parameter struct. Check against the desired FermAct name
			SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(const
					TwoFlavorMultihasenCancleMonomialParams& param_);

			// Copy Constructor
			SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial(const
					SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial& m):
				phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor){}

	//		//! even-even ln det Clover
			Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
					multi1d<LatticeColorMatrix> >& s){
				return Double(0);
			}

	  	//! odd-odd 
			Double S_odd_odd(const AbsFieldState<multi1d<LatticeColorMatrix>,
					multi1d<LatticeColorMatrix> >& s);

			void dsdq(P& F, const AbsFieldState<P,Q>& s);

			void refreshInternalFields(const AbsFieldState<P,Q>& field_state);

			void setInternalFields(const Monomial<P,Q>& m);

		protected:
			T& getPhi(void){
				return phi;
			}

			const T& getPhi(void) const{
				return phi;
			}
			// override base class pure virtual function
			virtual const SymEvenOddPrecWilsonTypeFermAct<T,P,Q>& getFermAct() const {};
			// Can't return const and can't use covariant return rule
			// so write a new one
			ShiftSymEvenOddPrecCloverFermAct& getFermAct(){
				return *fermact;
			}
		
			AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void){
				return *chrono_predictor;
			}

			const GroupXML_t& getInvParams(void) const {
				return inv_param;
			}

		private:
			// Hide empty constructor and =
			SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial();
			void operator=(const SymEvenOddPrecConstDetTwoFlavorWilsonMultihasenCancleMonomial&);

			// Pseudofermion field phi
			T phi;

			// A handle for the ShiftSymEvenOddPrecCloverFermAct
			Handle<ShiftSymEvenOddPrecCloverFermAct> fermact;
	
			// Shifted mass parameter
			Real mu;

			// The parameters for the inversion
			GroupXML_t inv_param;
			Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
	};

}

#endif
