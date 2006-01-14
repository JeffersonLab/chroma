// -*- C++ -*-
// $Id: prec_two_flavor_monomial5d_w.h,v 2.1 2006-01-14 05:22:32 edwards Exp $

/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#ifndef __prec_two_flavor_monomial5d_w_h__
#define __prec_two_flavor_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    extern const bool registered;
  };


  //! Wrapper class for 5D 2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
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
						   const TwoFlavorWilsonTypeFermMonomialParams& param_);

      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
//      EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(Handle< const EvenOddPrecWilsonFermAct >& fermact_, const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

#if 0
      const multi1d<LatticeFermion>& debugGetPhi(void) const {
	return getPhi();
      }


      void debugGetX(multi1d<LatticeFermion>& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s)  {
	getX(X,s);
      }
#endif

      const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
 
      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) {
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
      int getX(multi1d<LatticeFermion>& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s);


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
