// -*- C++ -*-
// $Id: prec_one_flavor_rat_monomial5d_w.h,v 1.1 2005-01-28 02:15:32 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial5d_w_h__
#define __prec_one_flavor_rat_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_w.h"

namespace Chroma 
{

  namespace EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  struct EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams {
    // Base Constructor
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams();

    // Read monomial from some root path
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
  };

  void read(XMLReader& xml, const string& path, EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams& params);


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D :
    public  OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial5D< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D(const string& fermact_name, 
						   const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams& param_);

      // Copy Constructor
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param) {}

      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) const {
	return Double(0);
      }


    protected:

      const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
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

      //! Return the numerator coefficient in force calc. partial fraction expansion
      const multi1d<Real>& getFPartFracCoeff() const {return FPartFracCoeff;}

      //! Return the denominator roots in force calc. partial fraction expansion
      const multi1d<Real>& getFPartFracRoot() const {return FPartFracRoot;}

      //! Return the numerator coefficient in force calc. partial fraction expansion for PV
      const multi1d<Real>& getFPVPartFracCoeff() const {return FPVPartFracCoeff;}

      //! Return the denominator roots in force calc. partial fraction expansion for PV
      const multi1d<Real>& getFPVPartFracRoot() const {return FPVPartFracRoot;}

      //! Return the numerator coefficient in heat-bath partial fraction expansion
      const multi1d<Real>& getHBPartFracCoeff() const {return HBPartFracCoeff;}

      //! Return the denominator roots in heat-bath partial fraction expansion
      const multi1d<Real>& getHBPartFracRoot() const {return HBPartFracRoot;}

      //! Return the numerator coefficient in heat-bath partial fraction expansion for PV
      const multi1d<Real>& getHBPVPartFracCoeff() const {return HBPVPartFracCoeff;}

      //! Return the denominator roots in heat-bath partial fraction expansion for PV
      const multi1d<Real>& getHBPVPartFracRoot() const {return HBPVPartFracRoot;}

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
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D();
      void operator=(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // Pseudofermion field phi
      multi1d<LatticeFermion> phiPV;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Coefficients and roots of partial fractions
      multi1d<Real> FPartFracCoeff;
      multi1d<Real> FPartFracRoot;
      multi1d<Real> FPVPartFracCoeff;
      multi1d<Real> FPVPartFracRoot;
      multi1d<Real> HBPartFracCoeff;
      multi1d<Real> HBPartFracRoot;
      multi1d<Real> HBPVPartFracCoeff;
      multi1d<Real> HBPVPartFracRoot;
    };


} //end namespace chroma

#endif
