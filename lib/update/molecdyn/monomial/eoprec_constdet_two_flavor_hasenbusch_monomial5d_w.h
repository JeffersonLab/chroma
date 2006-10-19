// -*- C++ -*-
// $Id: eoprec_constdet_two_flavor_hasenbusch_monomial5d_w.h,v 3.1 2006-10-19 16:01:34 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_hasenbusch_monomial5d_w_h__
#define __prec_two_flavor_hasenbusch_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D :
    public  TwoFlavorExactEvenOddPrecConstDetHasenbuschWilsonTypeFermMonomial5D< 
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
    EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D(const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param_);

    // Copy Constructor
    EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D(const EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact(m.fermact), fermact_prec(m.fermact_prec), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

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

    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getFermActPrec(void) const { 
      return *fermact_prec;
    }


    AbsChronologicalPredictor5D<T>& getMDSolutionPredictor(void) { 
      return *chrono_predictor;
    };

    //! Get parameters for the inverter
    const GroupXML_t getInvParams(void) const { 
      return inv_param;
    }

  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D();
    void operator=(const EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D&);

    // Pseudofermion field phi
    multi1d<T> phi;

    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact;

    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact_prec;

    // The parameters for the inversion
    GroupXML_t inv_param;

    Handle<AbsChronologicalPredictor5D<T> > chrono_predictor;
  };


} //end namespace chroma

#endif
