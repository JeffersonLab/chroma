// -*- C++ -*-
// $Id: prec_constdet_two_flavor_hasenbusch_monomial_w.h,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_hasenbusch_monomial_w_h__
#define __prec_two_flavor_hasenbusch_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomialEnv 
  {
    extern const std::string name;
    extern const bool registered;
  };


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecConstDetHasenbuschWilsonTypeFermMonomial< 
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
      EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial(const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param_);

      // Copy Constructor
      EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial(const EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), fermact_prec(m.fermact_prec), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

    protected:

      T& getPhi(void) {
	return phi;
      }

      const T& getPhi(void) const {
	return phi;
      }

      const EvenOddPrecWilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
	return *fermact;
      }

      const EvenOddPrecWilsonTypeFermAct<T,P,Q>& getFermActPrec(void) const { 
	return *fermact_prec;
      }


      AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void) { 
	return *chrono_predictor;
      };

    //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

    private:
 
      // Hide empty constructor and =
      EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial();
      void operator=(const EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      T phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact;

      Handle<const EvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact_prec;

      // The parameters for the inversion
      InvertParam_t inv_param;

      Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
    };


}; //end namespace chroma

#endif
