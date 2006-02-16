// -*- C++ -*-
// $Id: prec_constdet_one_flavor_rat_monomial5d_w.h,v 2.3 2006-02-16 02:59:03 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial5d_w_h__
#define __prec_one_flavor_rat_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    extern const std::string name;
    extern const bool registered;
  };


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
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(const OneFlavorWilsonTypeFermRatMonomial5DParams& param_);

      // Copy Constructor
      EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(const EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), nthRoot(m.nthRoot), nthRootPV(m.nthRootPV) {}

      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) {
	return Double(0);
      }


    protected:

      const EvenOddPrecConstDetWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const InvertParam_t getInvParams(void) const { 
	return inv_param;
      }

      //! Accessor for pseudofermion (read only)
      const multi1d< multi1d<LatticeFermion> >& getPhi(void) const {return phi;}

      //! mutator for pseudofermion
      multi1d< multi1d<LatticeFermion> >& getPhi(void) {return phi;}

      //! Accessor for PV pseudofermion (read only)
      const multi1d< multi1d<LatticeFermion> >& getPhiPV(void) const {return phiPV;}

      //! mutator for PV pseudofermion 
      multi1d< multi1d<LatticeFermion> >& getPhiPV(void) {return phiPV;}

      //! Return number of roots in used
      int getNthRoot() const {return nthRoot;}

      int getNthRootPV() const { return nthRootPV; }

      //! Return the partial fraction expansion for the action calc
      const RemezCoeff_t& getSPFE() const {return spfe;}

      //! Return the partial fraction expansion for the heat-bath
      const RemezCoeff_t& getSIPFE() const {return sipfe;}

      //! Return the partial fraction expansion for the force calc in PV
      const RemezCoeff_t& getFPVPFE() const {return fpvpfe;}

      //! Return the partial fraction expansion for the action calc in PV
      const RemezCoeff_t& getSPVPFE() const {return spvpfe;}

      //! Return the partial fraction expansion for the heat-bath in PV
      const RemezCoeff_t& getSIPVPFE() const {return sipvpfe;}

    private:
 
      // Hide empty constructor and =
      EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D();
      void operator=(const EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D&);

      // Pseudofermion field phi
      multi1d< multi1d<LatticeFermion> > phi;
      
      // Pseudofermion field phi
      multi1d< multi1d<LatticeFermion> > phiPV;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecConstDetWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Number of nth-roots
      int nthRoot;

      // Number of PV nth-roots
      int nthRootPV;

      //! Return the partial fraction expansion for the force calc
      const RemezCoeff_t& getFPFE() const {return fpfe;}
      // Coefficients and roots of partial fractions
      RemezCoeff_t  fpfe;
      RemezCoeff_t  spfe;
      RemezCoeff_t  sipfe;

      RemezCoeff_t  fpvpfe;
      RemezCoeff_t  spvpfe;
      RemezCoeff_t  sipvpfe;
    };


} //end namespace chroma

#endif
