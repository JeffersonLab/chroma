// -*- C++ -*-
// $Id: unprec_two_flavor_hasenbusch_monomial_w.h,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_two_flavor_hasenbusch_monomial_w_h__
#define __unprec_two_flavor_hasenbusch_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorHasenbuschWilsonTypeFermMonomialEnv 
  {
    extern const std::string name;
    extern const bool registered;
  };


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial 
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial :
    public  TwoFlavorExactUnprecHasenbuschWilsonTypeFermMonomial< 
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
      UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param_);


      // Copy Constructor
      UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(const UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial& m) : phi(m.phi), fermact((m.fermact)), fermact_prec(m.fermact_prec), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

    protected:

      T& getPhi(void) {
	// If phi are changed we must reset the chrono predictor
	return phi;
      }

      const T& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
	return *fermact;
      }

      const UnprecWilsonTypeFermAct<T,P,Q>& getFermActPrec(void) const { 
	return *fermact_prec;
      }

      AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void) {
	return *chrono_predictor;
      }

      //! Do an inversion of the type 
      const InvertParam_t getInvParams(void) const {
	return inv_param;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial();
      void operator=(const UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      T phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct<T,P,Q> > fermact;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct<T,P,Q> > fermact_prec;

      // The parameters for the inversion
      InvertParam_t inv_param;
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor4D<T> > chrono_predictor;
    };


} //end namespace chroma



#endif
