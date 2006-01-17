// -*- C++ -*-
// $Id: prec_constdet_one_flavor_rat_monomial_w.h,v 2.2 2006-01-17 16:01:46 bjoo Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial_w_h__
#define __prec_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    extern const bool registered;
  };


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial :
    public  OneFlavorRatExactEvenOddPrecConstDetWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial(const string& fermact_name, 
							    const OneFlavorWilsonTypeFermRatMonomialParams& param_);
      // Copy Constructor
      EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial(const EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial& m) 
	: phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), nthRoot(m.nthRoot) {}

    protected:

      multi1d<LatticeFermion>& getPhi(void) {return phi;}
      const multi1d<LatticeFermion>& getPhi(void) const {return phi;}

      const EvenOddPrecConstDetWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

      //! Return number of roots in used
      int getNthRoot() const {return nthRoot;}

      //! Return the partial fraction expansion for the force calc
      const RemezCoeff_t& getFPFE() const {return fpfe;}

      //! Return the partial fraction expansion for the action calc
      const RemezCoeff_t& getSPFE() const {return spfe;}

      //! Return the partial fraction expansion for the heat-bath
      const RemezCoeff_t& getSIPFE() const {return sipfe;}

    private:
 
      // Hide empty constructor and =
      EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial();
      void operator=(const EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecConstDetWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Number of nth-roots
      int nthRoot;

      // Coefficients and roots of partial fractions
      RemezCoeff_t  fpfe;
      RemezCoeff_t  spfe;
      RemezCoeff_t  sipfe;
    };


}; //end namespace chroma

#endif
