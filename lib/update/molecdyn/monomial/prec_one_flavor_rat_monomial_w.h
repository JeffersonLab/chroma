// -*- C++ -*-
// $Id: prec_one_flavor_rat_monomial_w.h,v 1.1 2005-01-28 02:15:32 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial_w_h__
#define __prec_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_w.h"

namespace Chroma 
{

  namespace EvenOddPrecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  struct EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams {
    // Base Constructor
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams();

    // Read monomial from some root path
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;

  };

  void read(XMLReader& xml, const string& path, EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& params);


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecOneFlavorWilsonTypeFermRatMonomial :
    public  OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(const string& fermact_name, 
						   const EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& param_);

      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
//      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(Handle< const EvenOddPrecWilsonFermAct >& fermact_, const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param) {}

      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) const {
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

      //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
      int getX(multi1d<LatticeFermion>& X, 
	       const multi1d<Real>& shifts, 
	       const LatticeFermion& chi, 
	       const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;

      //! Return the numerator coefficient in force calc. partial fraction expansion
      const multi1d<Real>& getFPartFracCoeff() const {return FPartFracCoeff;}

      //! Return the denominator roots in force calc. partial fraction expansion
      const multi1d<Real>& getFPartFracRoot() const {return FPartFracRoot;}

      //! Return the numerator coefficient in heat-bath partial fraction expansion
      const multi1d<Real>& getHBPartFracCoeff() const {return HBPartFracCoeff;}

      //! Return the denominator roots in heat-bath partial fraction expansion
      const multi1d<Real>& getHBPartFracRoot() const {return HBPartFracRoot;}

    private:
 
      // Hide empty constructor and =
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial();
      void operator=(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Coefficients and roots of partial fractions
      multi1d<Real> FPartFracCoeff;
      multi1d<Real> FPartFracRoot;
      multi1d<Real> HBPartFracCoeff;
      multi1d<Real> HBPartFracRoot;
    };


}; //end namespace chroma

#endif
