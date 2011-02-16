// -*- C++ -*-
// $Id: unprec_two_flavor_ratio_conv_conv_monomial_w.h,v 3.1 2008-05-23 21:31:36 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_two_flavor_ratio_conv_conv_monomial_w_h__
#define __unprec_two_flavor_ratio_conv_conv_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomialEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial 
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial :
    public  TwoFlavorExactUnprecRatioConvConvWilsonTypeFermMonomial< 
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
      UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial(const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param_);

    protected:

      T& getPhi() {
	// If phi are changed we must reset the chrono predictor
	return phi;
      }

      const T& getPhi() const {
	return phi;
      }

      const UnprecWilsonTypeFermAct<T,P,Q>& getNumerFermAct() const { 
	return *fermact_num;
      }

      const UnprecWilsonTypeFermAct<T,P,Q>& getDenomFermAct() const { 
	return *fermact_den;
      }

      AbsChronologicalPredictor4D<T>& getMDSolutionPredictor() {
	return *chrono_predictor;
      }

      //! Do an inversion of the type 
      const GroupXML_t& getNumerInvParams() const {
	return invParam_num;
      }

      //! Do an inversion of the type 
      const GroupXML_t& getDenomInvParams() const {
	return invParam_den;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial();
      void operator=(const UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      T phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct<T,P,Q> > fermact_num;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct<T,P,Q> > fermact_den;

      // The parameters for the inversion
      GroupXML_t invParam_num;
      GroupXML_t invParam_den;
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor4D<T> > chrono_predictor;
    };


} //end namespace chroma



#endif
