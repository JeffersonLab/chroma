// -*- C++ -*-
// $Id: prec_two_flavor_monomial_w.h,v 1.3 2005-04-10 21:46:42 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_monomial_w_h__
#define __prec_two_flavor_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv 
  {
    extern const bool registered;
  };

  // Parameter structure
  /*! @ingroup monomial */
  struct EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams {
    // Base Constructor
    EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
    string predictor_xml;   // The Chrono Predictor XML

  };

  void read(XMLReader& xml, const string& path, EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& params);


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a parameter struct. Check against the desired FermAct name
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial(const string& fermact_name, 
						   const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& param_);

      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
//      EvenOddPrecTwoFlavorWilsonTypeFermMonomial(Handle< const EvenOddPrecWilsonFermAct >& fermact_, const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial(const EvenOddPrecTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

#if 0
      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s)  {
	getX(X,s);
      }
#endif

      const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
 
      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) {
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

      // Do inversion M^dag M X = phi
      int getX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s);


      //! Get X = (A^dag*A)^{-1} eta
      int invert(LatticeFermion& X, 
		 const LinearOperator<LatticeFermion>& A,
		 const LatticeFermion& eta) const;

      AbsChronologicalPredictor4D<LatticeFermion>& getMDSolutionPredictor(void) { 
	return *chrono_predictor;
      };

    private:
 
      // Hide empty constructor and =
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial();
      void operator=(const EvenOddPrecTwoFlavorWilsonTypeFermMonomial&);

      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;

      Handle<AbsChronologicalPredictor4D<LatticeFermion> > chrono_predictor;
    };


}; //end namespace chroma

#endif
