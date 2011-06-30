// -*- C++ -*-
// $Id: lcm_integrator_leaps.h,v 3.4 2007-03-23 20:21:39 bjoo Exp $

#ifndef LCM_INTEGRATOR_LEAPS
#define LCM_INTEGRATOR_LEAPS

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "handle.h"

#include "singleton.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/integrator_shared.h"
namespace Chroma 
{

  //! LatticeColorMatrix integrator leaps
  /*! @ingroup integrator */
  namespace LCMMDIntegratorSteps 
  {
    class AnisoStepSizeArray {
    public:
      AnisoStepSizeArray() {
	step_factors.resize(Nd);
	for(int mu=0; mu < Nd; mu++) { 
	  step_factors[mu] = Real(1);
	}
      }
      
      inline void
      setAnisoStepSize(int t_dir_, const Real& step_factor_) {
	QDPIO::cout << "Setting dir " << t_dir_ << " to aniso factor " << step_factor_ << " originally " << Real(1)/(step_factor_*step_factor_) << endl;
	if( (t_dir_ >= 0) && (t_dir_ < Nd) ) {
	  step_factors[t_dir_] = step_factor_;
	}
	else {
	  QDPIO::cout << "Error t_dir must be between 0 and " << Nd-1 << ". It is " << t_dir_ << endl;
	  QDP_abort(1);
	}
      }

      inline Real 
      getStepSizeFactor(int t_dir) {
	if( (t_dir < 0) || (t_dir >= Nd) ) {
	  QDPIO::cout << "Error t_dir must be between 0 and " << Nd-1 << ". It is " << t_dir << endl;
	  QDP_abort(1);
	}
	return step_factors[t_dir];
      }

    private:
      multi1d<Real> step_factors;
    };

    typedef SingletonHolder<  AnisoStepSizeArray > theAnisoStepSizeArray;

    //! Leap with Q (with all monomials)
    /*! @ingroup integrator */
    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);


    // LeapP update for just a list of Monomials in an array of handles
    void leapP(const multi1d< IntegratorShared::MonomialPair >& monomials,
	                                       
	       const Real& dt, 
	       	       
	       AbsFieldState<multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >& s);




  } // End Namespace MDIntegratorSteps

} // End Namespace Chroma 

#endif
