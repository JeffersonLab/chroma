// -*- C++ -*-
// $Id: unprec_one_flavor_rat_monomial_w.h,v 1.4 2005-04-18 16:23:24 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_one_flavor_rat_monomial_w_h__
#define __unprec_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  /*! @ingroup monomial */
  struct UnprecOneFlavorWilsonTypeFermRatMonomialParams 
  {
    // Base Constructor
    UnprecOneFlavorWilsonTypeFermRatMonomialParams();

    // Read monomial from some root path
    UnprecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t   inv_param; // Inverter Parameters
    std::string     ferm_act;
    int             nthRoot;  // Use "n" copies of nth-root 1-flavor

    struct Remez_t   // eigenvalue bounds of M^dag*M
    {
      Real lowerMin;
      Real upperMax;
      int  forceDegree;
      int  actionDegree;
      int  digitPrecision;
    } remez;
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, UnprecOneFlavorWilsonTypeFermRatMonomialParams& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const UnprecOneFlavorWilsonTypeFermRatMonomialParams& params);

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
      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecOneFlavorWilsonTypeFermRatMonomial(const string& fermact_name, 
					      const UnprecOneFlavorWilsonTypeFermRatMonomialParams& param_);


      // Copy Constructor
      UnprecOneFlavorWilsonTypeFermRatMonomial(const UnprecOneFlavorWilsonTypeFermRatMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param), nthRoot(m.nthRoot) {}


    protected:

      multi1d<LatticeFermion>& getPhi(void) {return phi;}
      const multi1d<LatticeFermion>& getPhi(void) const {return phi;}

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
      int getX(multi1d<LatticeFermion>& X, 
	       const multi1d<Real>& shifts, 
	       const LatticeFermion& chi, 
	       const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;

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
      multi1d<LatticeFermion> phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Number of nth-roots
      int nthRoot;

      // Coefficients and roots of partial fractions
      RemezCoeff_t  fpfe;
      RemezCoeff_t  spfe;
      RemezCoeff_t  sipfe;
  };

} //end namespace chroma


#endif
