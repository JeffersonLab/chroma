// -*- C++ -*-
// $Id: unprec_two_flavor_wilson_monomial_w.h,v 1.1 2005-01-10 19:59:11 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef UNPREC_TWO_FLAVOR_WILSON_TYPE_MONOMIAL_H
#define UNPREC_TWO_FLAVOR_WILSON_TYPE_MONOMIAL_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"

namespace Chroma 
{

  namespace UnprecTwoFlavorWilsonTypeFermMonomialEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  struct UnprecTwoFlavorWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    UnprecTwoFlavorWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    UnprecTwoFlavorWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
  };

  void read(XMLReader& xml, const string& path, UnprecTwoFlavorWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const UnprecTwoFlavorWilsonTypeFermMonomialParams& params);

  //! Wrapper class for  2-flavor unprec ferm monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactUnprecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecTwoFlavorWilsonTypeFermMonomial(const string& fermact_name, 
					      const UnprecTwoFlavorWilsonTypeFermMonomialParams& param_);

      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
//      UnprecTwoFlavorWilsonTypeFermMonomial(const std::string& FermAct, 
//					      Handle< const UnprecWilsonTypeFermAct >& fermact_, 
//					      const InvertParam_t& inv_param_) : 
//	fermact(fermact_), inv_param(inv_param_) {init(FermAct);}

      // Copy Constructor
      UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param) {}

      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {
	getX(X,s);
      }

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      

    protected:

      LatticeFermion& getPhi(void) {
	return phi;
      }

      const LatticeFermion& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      //! Do inversion M^dag M X = phi
      void getX(LatticeFermion& X, 
		const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;

      //! Get X = (A^dag*A)^{-1} eta
      int invert(LatticeFermion& X, 
		 const LinearOperator<LatticeFermion>& A,
		 const LatticeFermion& eta) const;


    private:
      // Hide empty constructor and =
      UnprecTwoFlavorWilsonTypeFermMonomial();
      void operator=(const UnprecTwoFlavorWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
    };


} //end namespace chroma



#endif
