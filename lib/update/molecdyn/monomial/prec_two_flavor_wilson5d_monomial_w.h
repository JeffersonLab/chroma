// -*- C++ -*-
// $Id: prec_two_flavor_wilson5d_monomial_w.h,v 1.1 2005-01-13 15:55:04 bjoo Exp $

/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef PREC_TWO_FLAVOR_WILSON_TYPE_MONOMIAL5D_H
#define PREC_TWO_FLAVOR_WILSON_TYPE_MONOMIAL5D_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"

namespace Chroma 
{

  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  struct EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams {
    // Base Constructor
    EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams();

    // Read monomial from some root path
    EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
    string predictor_xml;
  };

  void read(XMLReader& xml, const string& path, EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams& params);


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D :
    public  TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(const string& fermact_name, 
						   const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams& param_);

      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
//      EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(Handle< const EvenOddPrecWilsonFermAct >& fermact_, const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

      const multi1d<LatticeFermion>& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(multi1d<LatticeFermion>& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s)  {
	getX(X,s);
      }

      const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
 
      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) const {
	return Double(0);
      }


    protected:

      multi1d<LatticeFermion>& getPhi(void) {
	return phi;
      }

      const multi1d<LatticeFermion>& getPhi(void) const {
	return phi;
      }

      const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      // Do inversion M^dag M X = phi
      int getX(multi1d<LatticeFermion>& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const ;


      //! Get X = (PV^dag*PV)^{-1} eta
      int getXPV(multi1d<LatticeFermion>& X, const multi1d<LatticeFermion>& eta, 
		  const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;


      //! Get X = (A^dag*A)^{-1} eta
      int invert(multi1d<LatticeFermion>& X, 
		 const LinearOperator< multi1d<LatticeFermion> >& A,
		 const multi1d<LatticeFermion>& eta) const;

      AbsChronologicalPredictor5D<LatticeFermion>& getMDSolutionPredictor(void) { 
	return *chrono_predictor;
      }

      
    private:
 
      // Hide empty constructor and =
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D();
      void operator=(const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D&);

      // Pseudofermion field phi
      multi1d<LatticeFermion> phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
      Handle<AbsChronologicalPredictor5D<LatticeFermion> > chrono_predictor;
    };


}; //end namespace chroma

#endif
