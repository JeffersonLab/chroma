// -*- C++ -*-
// $Id: unprec_one_flavor_rat_monomial_w.h,v 3.3 2008-01-23 18:23:36 bjoo Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_one_flavor_rat_monomial_w_h__
#define __unprec_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecOneFlavorWilsonTypeFermRatMonomial :
    public  OneFlavorRatExactUnprecWilsonTypeFermMonomial< 
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
      UnprecOneFlavorWilsonTypeFermRatMonomial(const OneFlavorWilsonTypeFermRatMonomialParams& param_);


      // Copy Constructor
      UnprecOneFlavorWilsonTypeFermRatMonomial(const UnprecOneFlavorWilsonTypeFermRatMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param), nthRoot(m.nthRoot) {}


    protected:

      multi1d<T>& getPhi(void) {return phi;}
      const multi1d<T>& getPhi(void) const {return phi;}

      const WilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const GroupXML_t& getInvParams(void) const { 
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
      UnprecOneFlavorWilsonTypeFermRatMonomial();
      void operator=(const UnprecOneFlavorWilsonTypeFermRatMonomial&);

      // Pseudofermion field phi
      multi1d<T> phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const WilsonTypeFermAct<T,P,Q> > fermact;

      // The parameters for the inversion
      GroupXML_t inv_param;

      // Number of nth-roots
      int nthRoot;

      // Coefficients and roots of partial fractions
      RemezCoeff_t  fpfe;
      RemezCoeff_t  spfe;
      RemezCoeff_t  sipfe;
  };

} //end namespace chroma


#endif
