// -*- C++ -*-
// $Id: unprec_two_flavor_monomial5d_w.h,v 2.2 2006-01-14 06:42:07 edwards Exp $

/*! @file
 * @brief Two-flavor collection of unpreconditioned 5D ferm monomials
 */

#ifndef __unprec_two_flavor_monomial5d_w_h__
#define __unprec_two_flavor_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    extern const bool registered;
  };


  //! Wrapper class for 5D 2-flavor unprec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorWilsonTypeFermMonomial5D :
    public  TwoFlavorExactUnprecWilsonTypeFermMonomial5D< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecTwoFlavorWilsonTypeFermMonomial5D(const string& fermact_name, 
					      const TwoFlavorWilsonTypeFermMonomialParams& param_);

      // Copy Constructor
      UnprecTwoFlavorWilsonTypeFermMonomial5D(const UnprecTwoFlavorWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

#if 0 
      const multi1d<LatticeFermion>& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(multi1d<LatticeFermion>& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {
	getX(X,s);
      }
#endif

      const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      

    protected:

      multi1d<LatticeFermion>& getPhi(void) {
	return phi;
      }

      const multi1d<LatticeFermion>& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }
     
      AbsChronologicalPredictor5D< LatticeFermion >& getMDSolutionPredictor(void) {
	return *chrono_predictor;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorWilsonTypeFermMonomial5D();
      void operator=(const UnprecTwoFlavorWilsonTypeFermMonomial5D&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Chrono Predictor
      Handle < AbsChronologicalPredictor5D<LatticeFermion> > chrono_predictor;
    };


} //end namespace chroma



#endif
