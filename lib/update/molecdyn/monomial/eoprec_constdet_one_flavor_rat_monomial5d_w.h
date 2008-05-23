// -*- C++ -*-
// $Id: eoprec_constdet_one_flavor_rat_monomial5d_w.h,v 3.2 2008-05-23 21:31:32 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial5d_w_h__
#define __prec_one_flavor_rat_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D :
    public  OneFlavorRatExactEvenOddPrecConstDetWilsonTypeFermMonomial5D< 
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
    EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(const OneFlavorWilsonTypeFermRatMonomialParams& param_);
    
    //! Even even contribution (eg ln det Clover)
    Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >& s) {
      return Double(0);
    }


  protected:

    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getFermAct(void) const { 
      return *fermact;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getActionInvParams(void) const { 
      return actionInvParam;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getForceInvParams(void) const { 
      return forceInvParam;
    }

    //! Accessor for pseudofermion (read only)
    const multi1d< multi1d<T> >& getPhi(void) const {return phi;}

    //! mutator for pseudofermion
    multi1d< multi1d<T> >& getPhi(void) {return phi;}

    //! Return number of roots in used
    int getNPF() const {return num_pf;}

    //! Return the partial fraction expansion for the action calc
    const RemezCoeff_t& getSPFE() const {return spfe;}

    //! Return the partial fraction expansion for the heat-bath
    const RemezCoeff_t& getSIPFE() const {return sipfe;}

  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D();
    void operator=(const EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D&);

    // Pseudofermion field phi
    multi1d< multi1d<T> > phi;
      
    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact;

    // The parameters for the inversion
    GroupXML_t actionInvParam;
    GroupXML_t forceInvParam;

    // Number of nth-roots
    int num_pf;

    //! Return the partial fraction expansion for the force calc
    const RemezCoeff_t& getFPFE() const {return fpfe;}
    // Coefficients and roots of partial fractions
    RemezCoeff_t  fpfe;
    RemezCoeff_t  spfe;
    RemezCoeff_t  sipfe;
  };


} //end namespace chroma

#endif
