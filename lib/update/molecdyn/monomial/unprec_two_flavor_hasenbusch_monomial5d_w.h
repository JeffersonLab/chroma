// -*- C++ -*-
// $Id: unprec_two_flavor_hasenbusch_monomial5d_w.h,v 3.3 2006-09-20 20:28:05 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_two_flavor_hasenbusch_monomial5d_w_h__
#define __unprec_two_flavor_hasenbusch_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial 
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D :
    public  TwoFlavorExactHasenbuschUnprecWilsonTypeFermMonomial5D< 
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix>,
      LatticeFermion> 
        {
    public: 
      // Typedefs to save typing
      typedef LatticeFermion      T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D(const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param_);

      // Copy Constructor
//      UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D(const UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact((m.fermact)), fermact_prec(m.fermact_prec), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

    protected:

      multi1d<T>& getPhi(void) {
	// If phi are changed we must reset the chrono predictor
	return phi;
      }

      const multi1d<T>& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct5D<T,P,Q>& getFermAct(void) const { 
	return *fermact;
      }

      const UnprecWilsonTypeFermAct5D<T,P,Q>& getFermActPrec(void) const { 
	return *fermact_prec;
      }

      AbsChronologicalPredictor5D<T>& getMDSolutionPredictor(void) {
	return *chrono_predictor;
      }

      //! Do an inversion of the type 
      const GroupXML_t getInvParams(void) const {
	return inv_param;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D();
      void operator=(const UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D&);

      // Pseudofermion field phi
      multi1d<T> phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact_prec;

      // The parameters for the inversion
      GroupXML_t inv_param;
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor5D<T> > chrono_predictor;
    };


} //end namespace chroma



#endif
