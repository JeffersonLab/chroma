// -*- C++ -*-
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __seoprec_two_flavor_ratio_conv_conv_monomial_w_h__
#define __seoprec_two_flavor_ratio_conv_conv_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace SymEvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomialEnv
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class SymEvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial :
    public  TwoFlavorExactSymEvenOddPrecConstDetRatioConvConvWilsonTypeFermMonomial<
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
    SymEvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial(const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param_);

  protected:

    T& getPhi() {
      return phi;
    }

    const T& getPhi() const {
      return phi;
    }

    const SymEvenOddPrecWilsonTypeFermAct<T,P,Q>& getNumerFermAct() const {
      return *fermact_num;
    }

    const SymEvenOddPrecWilsonTypeFermAct<T,P,Q>& getDenomFermAct() const {
      return *fermact_den;
    }


    AbsChronologicalPredictor4D<T>& getMDSolutionPredictor() { 
      return *chrono_predictor;
    };

    //! Get parameters for the inverter
    const GroupXML_t& getNumerInvParams() const { 
      return invParam_num;
    }

    const GroupXML_t& getDenomInvParams() const { 
      return invParam_den;
    }

  private:
 
    // Hide empty constructor and =
    SymEvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial();
    void operator=(const SymEvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial&);

    // Pseudofermion field phi
    T phi;

    // A handle for the SymEvenOddPrecWilsonFermAct
    Handle<const SymEvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact_num;

    Handle<const SymEvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact_den;

    // The parameters for the inversion
    GroupXML_t invParam_num;
    GroupXML_t invParam_den;

    Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
  };


} //end namespace chroma

#endif
