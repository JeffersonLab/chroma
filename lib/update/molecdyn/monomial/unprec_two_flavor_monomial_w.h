// -*- C++ -*-
// $Id: unprec_two_flavor_monomial_w.h,v 2.3 2006-01-14 05:56:59 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_two_flavor_monomial_w_h__
#define __unprec_two_flavor_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorWilsonTypeFermMonomialEnv 
  {
    extern const bool registered;
  };


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial 
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactUnprecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecTwoFlavorWilsonTypeFermMonomial(const string& fermact_name, 
					    const TwoFlavorWilsonTypeFermMonomialParams& param_);


      // Copy Constructor
      UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

#if 0
      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }


      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) {
	getX(X,s);
      }
#endif

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      

    protected:

      LatticeFermion& getPhi(void) {
	// If phi are changed we must reset the chrono predictor
	return phi;
      }

      const LatticeFermion& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

      AbsChronologicalPredictor4D<LatticeFermion>& getMDSolutionPredictor(void) {
	return *chrono_predictor;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorWilsonTypeFermMonomial();
      void operator=(const UnprecTwoFlavorWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor4D<LatticeFermion> > chrono_predictor;


    };


} //end namespace chroma



#endif
