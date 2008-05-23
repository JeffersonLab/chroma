// -*- C++ -*-
// $Id: unprec_two_flavor_monomial5d_w.h,v 3.3 2008-05-23 21:31:36 edwards Exp $

/*! @file
 * @brief Two-flavor collection of unpreconditioned 5D ferm monomials
 */

#ifndef __unprec_two_flavor_monomial5d_w_h__
#define __unprec_two_flavor_monomial5d_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace UnprecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    bool registerAll();
  }


  //! Wrapper class for 5D 2-flavor unprec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class UnprecTwoFlavorWilsonTypeFermMonomial5D :
    public  TwoFlavorExactUnprecWilsonTypeFermMonomial5D< 
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
      UnprecTwoFlavorWilsonTypeFermMonomial5D(const TwoFlavorWilsonTypeFermMonomialParams& param_);

    protected:

      multi1d<T>& getPhi(void) {
	return phi;
      }

      const multi1d<T>& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct5D<T,P,Q>& getFermAct(void) const { 
	return *fermact;
      }

      //! Get parameters for the inverter
      const GroupXML_t& getInvParams(void) const { 
	return inv_param;
      }
     
      AbsChronologicalPredictor5D<T>& getMDSolutionPredictor(void) {
	return *chrono_predictor;
      }

    private:
      // Hide empty constructor and =
      UnprecTwoFlavorWilsonTypeFermMonomial5D();
      void operator=(const UnprecTwoFlavorWilsonTypeFermMonomial5D&);

      // Pseudofermion field phi
      multi1d<T> phi;

      // A handle for the UnprecWilsonFermAct
      Handle<const UnprecWilsonTypeFermAct5D<T,P,Q> > fermact;

      // The parameters for the inversion
      GroupXML_t inv_param;

      // Chrono Predictor
      Handle < AbsChronologicalPredictor5D<T> > chrono_predictor;
    };


} //end namespace chroma



#endif
