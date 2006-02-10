// -*- C++ -*-
// $Id: prec_constdet_two_flavor_polynomial_monomial_w.h,v 2.2 2006-02-10 02:45:55 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_polynomial_monomial_w_h__
#define __prec_two_flavor_polynomial_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_polynomial_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv 
  {
    extern const bool registered;
    extern const std::string name;
  };


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecConstDetPolynomialWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial(const TwoFlavorWilsonTypeFermMonomialParams& param_);

      // Copy Constructor
      EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial(const EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param) {}

    protected:

      LatticeFermion& getPhi(void) {
	return phi;
      }

      const LatticeFermion& getPhi(void) const {
	return phi;
      }

      const PolyWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

    private:
 
      // Hide empty constructor and =
      EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial();
      void operator=(const EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const PolyWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
    };


}; //end namespace chroma

#endif
