// -*- C++ -*-
// $Id: prec_one_flavor_rat_monomial_w.h,v 1.6 2005-05-29 02:10:45 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_one_flavor_rat_monomial_w_h__
#define __prec_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/one_flavor_rat_monomial_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  /*! @ingroup monomial */
  struct EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams {
    // Base Constructor
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams();

    // Read monomial from some root path
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string& path);
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string& path,
							int expNumPower_, int expDenPower_);

    InvertParam_t   inv_param;     // Inverter Parameters
    std::string     ferm_act;
    int             expNumPower;   // (M^dag*M)^{expNumPower / (2*expDenPower)}
    int             expDenPower;   // (M^dag*M)^{expNumPower / (2*expDenPower)}
    int             nthRoot;       // Use "n" copies of nth-root 1-flavor

    struct Remez_t   // eigenvalue bounds of M^dag*M
    {
      Real lowerMin;
      Real upperMax;
      int  forceDegree;
      int  actionDegree;
      int  digitPrecision;
    } remez;
  };

  void read(XMLReader& xml, const string& path, EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& params);


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
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
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), nthRoot(m.nthRoot) {}

      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) {
	return Double(0);
      }


    protected:

      multi1d<LatticeFermion>& getPhi(void) {return phi;}
      const multi1d<LatticeFermion>& getPhi(void) const {return phi;}

      const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
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
      EvenOddPrecOneFlavorWilsonTypeFermRatMonomial();
      void operator=(const EvenOddPrecOneFlavorWilsonTypeFermRatMonomial&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      // Number of nth-roots
      int nthRoot;

      // Coefficients and roots of partial fractions
      RemezCoeff_t  fpfe;
      RemezCoeff_t  spfe;
      RemezCoeff_t  sipfe;
    };


}; //end namespace chroma

#endif
