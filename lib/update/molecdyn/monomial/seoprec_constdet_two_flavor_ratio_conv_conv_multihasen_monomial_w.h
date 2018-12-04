// -*- C++ -*-
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __SEOPREC_TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_W_H__
#define __SEOPREC_TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_W_H__
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_multihasen_monomial_params_w.h"
#include "actions/ferm/linop/shifted_seoprec_clover_linop_w.h"
#include "actions/ferm/fermacts/shifted_seoprec_clover_fermact_w.h"

namespace Chroma
{

	/*! @ingroup monomial */
	namespace SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialEnv
	{
		bool registerAll();
	}

	//! Wrapper class for  2-flavor even-odd prec ferm monomials
	/*! @ingroup monomial
	 *
	 * Monomial is expected to be the same for these fermacts
	 */
	class SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial:
		public ExactWilsonTypeFermMonomial<
		multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix>,
		LatticeFermion> 
	{
		public: 
			// Typedefs to save typing
			typedef LatticeFermion               T;
			typedef multi1d<LatticeColorMatrix>  P;
			typedef multi1d<LatticeColorMatrix>  Q;

			// Construct out of a parameter struct. Check against the desired FermAct name
			SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(const 
					TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param_);
			~SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(){}

			Double S(const AbsFieldState<P, Q>& s);

			// Total force of multi-hasen terms
			void dsdq(P& F, const AbsFieldState<P, Q>& s);

			void refreshInternalFields(const AbsFieldState<P, Q>& field_state);
			void setInternalFields(const Monomial<P, Q>& m);

			void resetPredictors(){
				getMDSolutionPredictor().reset();
			}

			int getNumHasenTerms() const {
				return numHasenTerms;
			}

		protected:

			// just override base class method which not needed at here
			//virtual const T& getPhi(void)const {};
			//virtual T& getPhi(void){};
			virtual const WilsonTypeFermAct<T,P,Q>& getFermAct(void)const{};
			// pesudo-fermion field for each Hasenbusch term
			T& getPhi(int i) {
				return phi[i];
			}
			const T& getPhi(int i) const {
				return phi[i];
			}

			ShiftSymEvenOddPrecCloverFermAct& getFermAct(){
				return *fermact;
			}

			AbsChronologicalPredictor4D<T>& getMDSolutionPredictor() { 
				return *chrono_predictor;
			};

			//! Get parameters for the inverter
			const GroupXML_t& getInvParams() const { 
				return invParam;
			}

		private:
			// Hide empty constructor and =
			SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial();
			void operator=(const SymEvenOddPrecConstDetTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial&);

			// Pseudofermion field phi for multi-hasenbusch term
			multi1d<T> phi;

			// A handle for the ShiftSymEvenOddPrecCloverFermAct
			Handle<ShiftSymEvenOddPrecCloverFermAct > fermact;

			// Shifted mass sym even-odd prec linear op
			//multi1d<Handle<ShiftSymEvenOddPrecCloverLinOp> > shift_linop;
			// Shifted mass sym even-odd prec ferm action
			//multi1d<Handle<const SymEvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact;
			// Shifted mass mu
			multi1d<Real> mu;
			// Number of Hasenbusch terms
			int numHasenTerms;

			// The parameters for the inversion
			GroupXML_t invParam;

			Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
	};


} //end namespace chroma
#endif
