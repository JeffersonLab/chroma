// -*- C++ -*-
// $Id: unprec_one_flavor_rat_monomial5d_w.h,v 1.2 2005-02-03 03:16:41 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 5D ferm monomials
 */

#ifndef __unprec_one_flavor_rat_monomial5d_w_h__
#define __unprec_one_flavor_rat_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_w.h"

namespace Chroma 
{

  namespace UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  struct UnprecOneFlavorWilsonTypeFermRatMonomial5DParams 
  {
    // Base Constructor
    UnprecOneFlavorWilsonTypeFermRatMonomial5DParams();

    // Read monomial from some root path
    UnprecOneFlavorWilsonTypeFermRatMonomial5DParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;

    struct Remez_t   // eigenvalue bounds of M^dag*M
    {
      Real lowerMin;
      Real upperMax;
      Real lowerMinPV;
      Real upperMaxPV;
      int  forceDegree;
      int  actionDegree;
      int  digitPrecision;
    } remez;
  };

  void read(XMLReader& xml, const string& path, UnprecOneFlavorWilsonTypeFermRatMonomial5DParams& param);

  void write(XMLWriter& xml, const string& path, const UnprecOneFlavorWilsonTypeFermRatMonomial5DParams& params);

  //! Wrapper class for 5D 2-flavor unprec ferm monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecOneFlavorWilsonTypeFermRatMonomial5D :
    public  OneFlavorRatExactUnprecWilsonTypeFermMonomial5D< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecOneFlavorWilsonTypeFermRatMonomial5D(const string& fermact_name, 
						 const UnprecOneFlavorWilsonTypeFermRatMonomial5DParams& param_);

      // Copy Constructor
      UnprecOneFlavorWilsonTypeFermRatMonomial5D(const UnprecOneFlavorWilsonTypeFermRatMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param) {}


    protected:

      const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Accessor for pseudofermion (read only)
      const multi1d<LatticeFermion>& getPhi(void) const {return phi;}

      //! mutator for pseudofermion
      multi1d<LatticeFermion>& getPhi(void) {return phi;}

      //! Accessor for PV pseudofermion (read only)
      const multi1d<LatticeFermion>& getPhiPV(void) const {return phiPV;}

      //! mutator for PV pseudofermion 
      multi1d<LatticeFermion>& getPhiPV(void) {return phiPV;}

      //! Return the partial fraction expansion for the force calc
      const RemezCoeff_t& getFPFE() const {return fpfe;}

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

      //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
      int getX(multi1d< multi1d<LatticeFermion> >& X, 
	       const multi1d<Real>& shifts, 
	       const multi1d<LatticeFermion>& chi, 
	       const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;

      //! Multi-mass solver  (V^dagV + q_i)^{-1} chi for PV using partfrac
      int getXPV(multi1d< multi1d<LatticeFermion> >& X, 
		 const multi1d<Real>& shifts, 
		 const multi1d<LatticeFermion>& chi, 
		 const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;
      
      //! Get X = (A^dag*A + q_i)^{-1} eta
      int invert(multi1d< multi1d<LatticeFermion> >& X, 
		 const multi1d<Real>& shifts, 
		 const LinearOperator< multi1d<LatticeFermion> >& A,
		 const multi1d<LatticeFermion>& eta) const;
     
    private:
      // Hide empty constructor and =
      UnprecOneFlavorWilsonTypeFermRatMonomial5D();
      void operator=(const UnprecOneFlavorWilsonTypeFermRatMonomial5D&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // Pseudofermion field phi
      multi1d<LatticeFermion> phiPV;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

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
