// -*- C++ -*-
// $Id: unprec_one_flavor_rat_monomial_w.h,v 1.1 2005-01-28 02:15:33 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_one_flavor_rat_monomial_w_h__
#define __unprec_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_w.h"

namespace Chroma 
{

  namespace UnprecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  struct UnprecOneFlavorWilsonTypeFermRatMonomialParams 
  {
    // Base Constructor
    UnprecOneFlavorWilsonTypeFermRatMonomialParams();

    // Read monomial from some root path
    UnprecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
  };

  void read(XMLReader& xml, const string& path, UnprecOneFlavorWilsonTypeFermRatMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const UnprecOneFlavorWilsonTypeFermRatMonomialParams& params);

  //! Wrapper class for  2-flavor unprec ferm monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecOneFlavorWilsonTypeFermRatMonomial :
    public  OneFlavorRatExactUnprecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecOneFlavorWilsonTypeFermRatMonomial(const string& fermact_name, 
					      const UnprecOneFlavorWilsonTypeFermRatMonomialParams& param_);


      // Copy Constructor
      UnprecOneFlavorWilsonTypeFermRatMonomial(const UnprecOneFlavorWilsonTypeFermRatMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param) {}


    protected:

      LatticeFermion& getPhi(void) {return phi;}

      const LatticeFermion& getPhi(void) const {return phi;}

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
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
      UnprecOneFlavorWilsonTypeFermRatMonomial();
      void operator=(const UnprecOneFlavorWilsonTypeFermRatMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Coefficients and roots of partial fractions
      multi1d<Real> FPartFracCoeff;
      multi1d<Real> FPartFracRoot;
      multi1d<Real> HBPartFracCoeff;
      multi1d<Real> HBPartFracRoot;
  };


} //end namespace chroma



#endif
