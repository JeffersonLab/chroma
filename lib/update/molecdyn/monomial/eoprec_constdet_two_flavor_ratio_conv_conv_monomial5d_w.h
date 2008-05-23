// -*- C++ -*-
// $Id: eoprec_constdet_two_flavor_ratio_conv_conv_monomial5d_w.h,v 3.1 2008-05-23 21:31:33 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_ratio_conv_conv_monomial5d_w_h__
#define __prec_two_flavor_ratio_conv_conv_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_ratio_conv_conv_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial5DEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial5D :
    public  TwoFlavorExactEvenOddPrecConstDetRatioConvConvWilsonTypeFermMonomial5D< 
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
    EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial5D(const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param_);

  protected:

    multi1d<T>& getPhi() {
      return phi;
    }

    const multi1d<T>& getPhi() const {
      return phi;
    }

    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getNumerFermAct() const { 
      return *fermact_num;
    }

    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getDenomFermAct() const { 
      return *fermact_den;
    }


    AbsChronologicalPredictor5D<T>& getMDSolutionPredictor() { 
      return *chrono_predictor;
    };

    //! Get parameters for the inverter
    const GroupXML_t getNumerInvParams() const { 
      return invParam_num;
    }

  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial5D();
    void operator=(const EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial5D&);

    // Pseudofermion field phi
    multi1d<T> phi;

    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact_num;

    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact_den;

    // The parameters for the inversion
    GroupXML_t invParam_num;

    Handle<AbsChronologicalPredictor5D<T> > chrono_predictor;
  };


} //end namespace chroma

#endif
