// -*- C++ -*-
// $Id: prec_constdet_two_flavor_monomial_w.h,v 2.4 2006-02-25 22:48:58 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_monomial_w_h__
#define __prec_two_flavor_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomialEnv 
  {
    extern const std::string name;
    extern const bool registered;
  };


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecConstDetWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial(const TwoFlavorWilsonTypeFermMonomialParams& param_);

      // Copy Constructor
      EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}


      const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
 
      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) {
	return Double(0);
      }


    protected:

      LatticeFermion& getPhi(void) {
	return phi;
      }

      const LatticeFermion& getPhi(void) const {
	return phi;
      }

      const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }
      
      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

      AbsChronologicalPredictor4D<LatticeFermion>& getMDSolutionPredictor(void) { 
	return *chrono_predictor;
      };


    private:
 
      // Hide empty constructor and =
      EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial();
      void operator=(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      Handle<AbsChronologicalPredictor4D<LatticeFermion> > chrono_predictor;
    };


}; //end namespace chroma

#endif
