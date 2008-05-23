// -*- C++ -*-
// $Id: eoprec_constdet_two_flavor_polyprec_monomial_w.h,v 3.2 2008-05-23 21:31:33 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_polyprec_monomial_w_h__
#define __prec_two_flavor_polyprec_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_polyprec_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomialEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecConstDetPolyPrecWilsonTypeFermMonomial< 
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
    EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial(const TwoFlavorWilsonTypeFermMonomialParams& param_);

    // Copy Constructor
    EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial(const EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

  protected:

    T& getPhi(void) {
      return phi;
    }

    const T& getPhi(void) const {
      return phi;
    }

    const PolyWilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
      return *fermact;
    }

    AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void) { 
      return *chrono_predictor;
    };

    //! Get parameters for the inverter
    const GroupXML_t& getInvParams(void) const { 
      return inv_param;
    }

  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial();
    void operator=(const EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial&);

    // Pseudofermion field phi
    T phi;

    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const PolyWilsonTypeFermAct<T,P,Q> > fermact;

    // The parameters for the inversion
    GroupXML_t inv_param;

    Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
  };


} //end namespace chroma

#endif
