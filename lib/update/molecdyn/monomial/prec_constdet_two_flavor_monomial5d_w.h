// -*- C++ -*-
// $Id: prec_constdet_two_flavor_monomial5d_w.h,v 2.2 2006-01-17 16:01:46 bjoo Exp $

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
    extern const bool registered;
  };


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
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D(const string& fermact_name, 
						   const TwoFlavorWilsonTypeFermMonomialParams& param_);

      // Copy Constructor
      EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

      const EvenOddPrecConstDetWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
    protected:

      multi1d<LatticeFermion>& getPhi(void) {
	return phi;
      }

      const multi1d<LatticeFermion>& getPhi(void) const {
	return phi;
      }

      const EvenOddPrecConstDetWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

      AbsChronologicalPredictor5D<LatticeFermion>& getMDSolutionPredictor(void) { 
	return *chrono_predictor;
      }

      
    private:
 
      // Hide empty constructor and =
      EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D();
      void operator=(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
      Handle<AbsChronologicalPredictor5D<LatticeFermion> > chrono_predictor;
    };


}; //end namespace chroma

#endif
