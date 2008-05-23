// -*- C++ -*-
// $Id: eoprec_constdet_two_flavor_ratio_conv_rat_monomial_w.h,v 3.1 2008-05-23 21:31:33 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __eoprec_two_flavor_ratio_conv_rat_monomial_w_h__
#define __eoprec_two_flavor_ratio_conv_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_rat_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_rat_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorRatioConvRatWilsonTypeFermMonomialEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorRatioConvRatWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecConstDetRatioConvRatWilsonTypeFermMonomial< 
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
    EvenOddPrecConstDetTwoFlavorRatioConvRatWilsonTypeFermMonomial(const TwoFlavorRatioConvRatWilsonTypeFermMonomialParams& param_);

  protected:

    T& getPhi() {
      return phi;
    }

    const T& getPhi() const {
      return phi;
    }

    const EvenOddPrecWilsonTypeFermAct<T,P,Q>& getNumerFermAct() const { 
      return *fermact_num;
    }

    const EvenOddPrecWilsonTypeFermAct<T,P,Q>& getDenomFermAct() const { 
      return *fermact_den;
    }


    AbsChronologicalPredictor4D<T>& getMDSolutionPredictor() { 
      return *chrono_predictor;
    };

    //! Get parameters for the inverter
    const GroupXML_t& getNumerInvParams() const { 
      return invParam_num;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getDenomActionInvParams(void) const { 
      return actionInvParam_den;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getDenomForceInvParams(void) const { 
      return forceInvParam_den;
    }

    //! Return the partial fraction expansion for the force calc
    const RemezCoeff_t& getDenomFPFE() const {return fpfe_den;}

    //! Return the partial fraction expansion for the action calc
    const RemezCoeff_t& getDenomSPFE() const {return spfe_den;}

    //! Return the partial fraction expansion for the heat-bath
    const RemezCoeff_t& getDenomSIPFE() const {return sipfe_den;}

  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetTwoFlavorRatioConvRatWilsonTypeFermMonomial();
    void operator=(const EvenOddPrecConstDetTwoFlavorRatioConvRatWilsonTypeFermMonomial&);

    // Pseudofermion field phi
    T phi;

    Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;

    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const EvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact_num;

    Handle<const EvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact_den;

    // The parameters for the inversion
    GroupXML_t invParam_num;

    // The parameters for the inversion
    GroupXML_t actionInvParam_den;
    GroupXML_t forceInvParam_den;

    // Coefficients and roots of partial fractions
    RemezCoeff_t  fpfe_den;
    RemezCoeff_t  spfe_den;
    RemezCoeff_t  sipfe_den;
  };


} //end namespace chroma

#endif
