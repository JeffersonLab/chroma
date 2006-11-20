// -*- C++ -*-
// $Id: lcm_integrator_leaps.h,v 3.2 2006-11-20 19:15:02 bjoo Exp $

#ifndef LCM_INTEGRATOR_LEAPS
#define LCM_INTEGRATOR_LEAPS

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "handle.h"

#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"

namespace Chroma 
{

  //! LatticeColorMatrix integrator leaps
  /*! @ingroup integrator */
  namespace LCMMDIntegratorSteps 
  {

    //! LeapP for just a selected list of monomials
    /*! @ingroup integrator */
    void leapP(const Real& dt, 

	       AbsHamiltonian<multi1d<LatticeColorMatrix>,
	                      multi1d<LatticeColorMatrix> >& H,

	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);

    //! Leap with Q (with all monomials)
    /*! @ingroup integrator */
    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);


    //! LeapP for just a selected list of monomials
    /*! @ingroup integrator */
    void leapP(const multi1d<int>& monomial_list,
	       const Real& dt, 

	        AbsHamiltonian<multi1d<LatticeColorMatrix>,
	                      multi1d<LatticeColorMatrix> >& H,

	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);

    // LeapP update for just a list of Monomials in an array of handles
    void leapP(const multi1d< 
	       Handle<Monomial< multi1d<LatticeColorMatrix>, 
	                        multi1d<LatticeColorMatrix> > >
	       > monomials,
	                                       
	       const Real& dt, 
	       	       
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s); 


  } // End Namespace MDIntegratorSteps

} // End Namespace Chroma 

#endif
