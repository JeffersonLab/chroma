// -*- C++ -*-
// $Id: unprec_two_flavor_wilson_monomial_w.h,v 1.2 2005-01-11 15:28:03 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef UNPREC_TWO_FLAVOR_WILSON_TYPE_MONOMIAL_H
#define UNPREC_TWO_FLAVOR_WILSON_TYPE_MONOMIAL_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "update/molecdyn/chrono_predictor.h"

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
    std::string predictor_xml;  // ChronologicalPredictor XML
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


      // Copy Constructor
      UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const  {
	getX(X,s);
      }

      const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      

    protected:

      LatticeFermion& getPhi(void) {
	// If phi are changed we must reset the chrono predictor
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
		const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const ;

      //! Get X = (A^dag*A)^{-1} eta
      int invert(LatticeFermion& X, 
		 const LinearOperator<LatticeFermion>& A,
		 const LatticeFermion& eta) const;

      AbsChronologicalPredictor4D<LatticeFermion>& getMDSolutionPredictor(void) {
	return *chrono_predictor;
      }

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
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor4D<LatticeFermion> > chrono_predictor;


    };


} //end namespace chroma



#endif
