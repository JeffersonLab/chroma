// -*- C++ -*-
// $Id: eoprec_constdet_one_flavor_ratio_rat_conv_monomial5d_w.h,v 3.1 2008-05-23 21:31:32 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_one_flavor_ratio_rat_conv_monomial5d_w_h__
#define __prec_one_flavor_ratio_rat_conv_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_ratio_rat_conv_monomial5d_w.h"
#include "update/molecdyn/monomial/one_flavor_ratio_rat_conv_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatRatMonomial5DEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatConvMonomial5D :
    public  OneFlavorRatioRatConvExactEvenOddPrecConstDetWilsonTypeFermMonomial5D< 
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
    EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatConvMonomial5D(const OneFlavorWilsonTypeFermRatioRatConvMonomialParams& param_);
    
    //! Even even contribution (eg ln det Clover)
    Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >& s) {
      return Double(0);
    }


  protected:

    //! Get at fermion action
    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getNumerFermAct() const { 
      return *fermact_num;
    }

    //! Get at fermion action
    const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>& getDenomFermAct() const { 
      return *fermact_den;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getNumerActionInvParams() const { 
      return actionInvParam_num;
    }

    //! Get parameters for the inverter
    const GroupXML_t& getNumerForceInvParams() const { 
      return forceInvParam_num;
    }

    //! Return number of roots in used
    int getNPF() const {return num_pf;}

    //! Accessor for pseudofermion (read only)
    const multi1d< multi1d<T> >& getPhi() const {return phi;}

    //! mutator for pseudofermion
    multi1d< multi1d<T> >& getPhi() {return phi;}

    //! Return the partial fraction expansion for the force calc
    const RemezCoeff_t& getNumerFPFE() const {return fpfe_num;}

    //! Return the partial fraction expansion for the action calc
    const RemezCoeff_t& getNumerSPFE() const {return spfe_num;}

    //! Return the partial fraction expansion for the heat-bath
    const RemezCoeff_t& getNumerSIPFE() const {return sipfe_num;}

  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatConvMonomial5D();
    void operator=(const EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatConvMonomial5D&);

    // Pseudofermion field phi
    multi1d< multi1d<T> > phi;
      
    // A handle for the EvenOddPrecConstDetWilsonFermAct
    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact_num;
    Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q> > fermact_den;

    // The parameters for the inversion
    GroupXML_t actionInvParam_num;
    GroupXML_t forceInvParam_num;

    // Number of nth-roots
    int num_pf;

    // Coefficients and roots of partial fractions
    RemezCoeff_t  fpfe_num;
    RemezCoeff_t  spfe_num;
    RemezCoeff_t  sipfe_num;
  };


} //end namespace chroma

#endif
