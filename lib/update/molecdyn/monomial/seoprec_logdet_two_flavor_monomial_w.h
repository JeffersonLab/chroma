// -*- C++ -*-
/*! @file
 * @brief Two-flavor collection of symmetric even-odd preconditioned 4D ferm monomials
 */

#ifndef __seoprec_two_flavor_logdet_monomial_w_h__
#define __seoprec_two_flavor_logdet_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomialEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor symmetric even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactSymEvenOddPrecLogDetWilsonTypeFermMonomial< 
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
    SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial(const TwoFlavorWilsonTypeFermMonomialParams& param_);

    // Construct from a fermact handle and inv params
    // FermAct already holds BC-s
//      SymEvenOddPrecTwoFlavorWilsonTypeFermMonomial(Handle< const SymEvenOddPrecWilsonFermAct >& fermact_, const GroupXML_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

    // Copy Constructor
    SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial(const SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

  protected:

    T& getPhi(void) {
      return phi;
    }

    const T& getPhi(void) const {
      return phi;
    }

    const SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
      return *fermact;
    }
      
    //! Get parameters for the inverter
    const GroupXML_t& getInvParams(void) const { 
      return inv_param;
    }

    AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void) { 
      return *chrono_predictor;
    };


  private:
 
    // Hide empty constructor and =
    SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial();
    void operator=(const SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial&);

    // Pseudofermion field phi
    T phi;

    // A handle for the SymEvenOddPrecWilsonFermAct
    Handle<const SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q> > fermact;

    // The parameters for the inversion
    GroupXML_t inv_param;

    Handle< AbsChronologicalPredictor4D<T> > chrono_predictor;
  };


} //end namespace chroma

#endif
