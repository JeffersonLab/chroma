
// -*- C++ -*-

/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef _UNPREC_TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_W_H__
#define _UNPREC_TWO_FLAVOR_RATIO_CONV_CONV_NULTIHASEN_MONOMIAL_W_H__

#include "update/molecdyn/monomial/unprec_two_flavor_ratio_conv_conv_monomial_w.h"

namespace Chroma
{
	namespace UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialEnv
	{
		bool registerAll();
	}

	class UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial
	{
		public:

			// Typedefs to save typing
			typedef LatticeFermion               T;
			typedef multi1d<LatticeColorMatrix>  P;
			typedef multi1d<LatticeColorMatrix>  Q;
			
			UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(const
					TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param_);
			~UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(){}

			Double S(const AbsFieldState<P,Q>& s);

			void dsdq(P& F, const AbsFieldState<P,Q>& s);

			void refreshInternalFields(const AbsFieldState<P,Q>& field_state);
			void setInternalFields(const Monomial<P,Q>& m);

			void resetPredictors(){
				getMDSolutionPredictor().reset();
			}
			int getNumHasenTerms(){
				return numHasenTerms;
			}
		private:
			int numHasenTerms;
			mult1d<Handle<UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial> > multihasenMonomial;

	}

}

#endif
