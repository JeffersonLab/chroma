// -*- C++ -*-
// $Id: prec_one_flavor_rat_monomial5d_w.h,v 1.6 2005-05-02 09:35:58 bjoo Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial5d_w_h__
#define __prec_one_flavor_rat_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  /*! @ingroup monomial */
  struct EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams 
  {
    // Base Constructor
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams();

    // Read monomial from some root path
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams(XMLReader& in, const std::string&  path);
    InvertParam_t   inv_param; // Inverter Parameters
    std::string     ferm_act;
    int             nthRoot;  // Use "n" copies of nth-root 1-flavor
    int             nthRootPV; // Use "n" copies of nth-root 1 flavor PV

    struct Remez_t   // eigenvalue bounds of M^dag*M
    {
      Real lowerMin;
      Real upperMax;
      Real lowerMinPV;
      Real upperMaxPV;

      int  degree;
      int  degreePV;
      
      int digitPrecision;
    } remez;
  };

  void read(XMLReader& xml, const string& path, EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DParams& params);


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
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
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), nthRoot(m.nthRoot), nthRootPV(m.nthRootPV) {}

      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) {
	return Double(0);
      }


    protected:

      const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
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
      multi1d< multi1d<LatticeFermion> > phi;
      
      // Pseudofermion field phi
      multi1d< multi1d<LatticeFermion> > phiPV;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

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
