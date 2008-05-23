// -*- C++ -*-
// $Id: eoprec_constdet_two_flavor_monomial5d_w.h,v 3.2 2008-05-23 21:31:33 edwards Exp $

/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_two_flavor_monomial5d_w_h__
#define __prec_two_flavor_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D :
    public  TwoFlavorExactEvenOddPrecConstDetWilsonTypeFermMonomial5D< 
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
    EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D(const TwoFlavorWilsonTypeFermMonomialParams& param_);

    // Copy Constructor
    EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

  protected:

    multi1d<T>& getPhi(void) {
      return phi;
    }

    const multi1d<T>& getPhi(void) const {
      return phi;
    }

    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getFermAct(void) const { 
      return *fermact;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getInvParams(void) const { 
      return inv_param;
    }

    AbsChronologicalPredictor5D<T>& getMDSolutionPredictor(void) { 
      return *chrono_predictor;
    }

      
  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D();
    void operator=(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D&);

    // Pseudofermion field phi
    multi1d<T> phi;

    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact;

    // The parameters for the inversion
    GroupXML_t inv_param;
    Handle<AbsChronologicalPredictor5D<T> > chrono_predictor;
  };


} //end namespace chroma

#endif
