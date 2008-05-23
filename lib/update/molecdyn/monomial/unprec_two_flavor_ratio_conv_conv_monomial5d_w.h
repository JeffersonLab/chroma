// -*- C++ -*-
// $Id: unprec_two_flavor_ratio_conv_conv_monomial5d_w.h,v 3.1 2008-05-23 21:31:36 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_two_flavor_ratio_conv_conv_monomial5d_w_h__
#define __unprec_two_flavor_ratio_conv_conv_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5DEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial 
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D :
    public  TwoFlavorExactUnprecRatioConvConvWilsonTypeFermMonomial5D< 
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix>,
      LatticeFermion> 
        {
    public: 
      // Typedefs to save typing
      typedef LatticeFermion      T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D(const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param_);

    protected:

      multi1d<T>& getPhi() {
	// If phi are changed we must reset the chrono predictor
	return phi;
      }

      const multi1d<T>& getPhi() const {
	return phi;
      }

      const UnprecWilsonTypeFermAct5D<T,P,Q>& getNumerFermAct() const { 
	return *fermact_num;
      }

      const UnprecWilsonTypeFermAct5D<T,P,Q>& getDenomFermAct() const { 
	return *fermact_den;
      }

      AbsChronologicalPredictor5D<T>& getMDSolutionPredictor() {
	return *chrono_predictor;
      }

      //! Do an inversion of the type 
      const GroupXML_t getNumerInvParams() const {
	return invParam_num;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D();
      void operator=(const UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D&);

      // Pseudofermion field phi
      multi1d<T> phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact_num;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact_den;

      // The parameters for the inversion
      GroupXML_t invParam_num;
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor5D<T> > chrono_predictor;
    };


} //end namespace chroma



#endif
