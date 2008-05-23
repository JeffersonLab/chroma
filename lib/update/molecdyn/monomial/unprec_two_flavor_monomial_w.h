// -*- C++ -*-
// $Id: unprec_two_flavor_monomial_w.h,v 3.4 2008-05-23 21:31:36 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#ifndef __unprec_two_flavor_monomial_w_h__
#define __unprec_two_flavor_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorWilsonTypeFermMonomialEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for  2-flavor unprec ferm monomials
  /*! @ingroup monomial 
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactUnprecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Typedefs to save typing
      typedef LatticeFermion               T;
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      // Construct out of a parameter struct. Check against the desired FermAct name
      UnprecTwoFlavorWilsonTypeFermMonomial(const TwoFlavorWilsonTypeFermMonomialParams& param_);


      // Copy Constructor
      UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact((m.fermact)), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}

    protected:

      T& getPhi(void) {
	// If phi are changed we must reset the chrono predictor
	return phi;
      }

      const T& getPhi(void) const {
	return phi;
      }

      const WilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const GroupXML_t& getInvParams(void) const { 
	return inv_param;
      }

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
      Handle<const WilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      GroupXML_t inv_param;
      
      // A handle for the chrono predictor
      Handle< AbsChronologicalPredictor4D<LatticeFermion> > chrono_predictor;


    };


} //end namespace chroma



#endif
