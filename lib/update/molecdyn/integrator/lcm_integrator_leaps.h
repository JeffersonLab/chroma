// -*- C++ -*-
// $Id: lcm_integrator_leaps.h,v 3.0 2006-04-03 04:59:07 edwards Exp $

#ifndef LCM_INTEGRATOR_LEAPS
#define LCM_INTEGRATOR_LEAPS

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"

namespace Chroma 
{

  //! LatticeColorMatrix integrator leaps
  /*! @ingroup integrator */
  namespace LCMMDIntegratorSteps 
  {

    //! LeapP for just a selected list of monomials
    void leapP(const Real& dt, 

	       AbsHamiltonian<multi1d<LatticeColorMatrix>,
	                      multi1d<LatticeColorMatrix> >& H,

	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);

    //! Leap with Q (with all monomials)
    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);


    //! LeapP for just a selected list of monomials
    void leapP(const multi1d<int>& monomial_list,
	       const Real& dt, 

	        AbsHamiltonian<multi1d<LatticeColorMatrix>,
	                      multi1d<LatticeColorMatrix> >& H,

	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);


  }; // End Namespace MDIntegratorSteps




}; // End Namespace Chroma 

#endif
