#ifndef ABS_INTEGRATOR_H
#define ABS_INTEGRATOR_H

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma {

  //! MD integrator interface
  template<typename P, typename Q>
  class  AbsMDIntegrator {
  public:

    // Virtual destructor
    virtual ~AbsMDIntegrator(void) {}

    // Do a trajectory -- that is all HMD/HMC needs to know
    virtual void operator()(AbsFieldState<P,Q>& s) = 0;

    //! Get at the MD Hamiltonian
    //  This needs to be exposed to initialise the PF fields.
    virtual AbsHamiltonian<P,Q>& getHamiltonian(void)  = 0;

  }; // End class MD integrator


  //! MD integrator interface for PQP leapfrog
  template<typename P, typename Q>
  class PQPLeapfrogIntegrator : public AbsMDIntegrator<P,Q> {
  public:
    //! Virtual destructor
    virtual ~PQPLeapfrogIntegrator(void) {} 


    //! Perform a trajectory 
    virtual void operator()(AbsFieldState<P,Q>& s) {
      Real dt = getStepSize(); 
      Real dtby2 = dt / Real(2);
      Real tau0 = getTrajLength();

      // Start time = 0
      Real t = Real(0);

      bool endP = false;

      // First half step by leapP
      leapP(dtby2, s);
	
      while(! endP ) { 

	leapQ(dt, s);
	
	t += dt;
	
	// Check if this was the last gauge update
	// Tricky because you cannot really check for t = tau
	// with any certainty due to potential FP equality testing issues
	// However if the time remaining is substantially less than dt
	// ie dtby2 is less than dt and is less than dt by any reasonable
	// epsilon then we finish
	if( toBool( fabs(tau0 - t) <  dtby2  ) ) {
	  // Time left is less than dtby2
	  // Finish with a half P leap and signal end
	  leapP(dtby2, s);
	  endP = true;
	  
	}
	else {
	  leapP(dt, s);
	}
      } // end while
    } // end function

    //! Get at the MD Hamiltonian
    //  This needs to be exposed to initialise the PF fields.
    virtual AbsHamiltonian<P,Q>& getHamiltonian(void)  = 0;
  
  protected:
    //! Leap with P
    virtual void leapP(const Real& dt, AbsFieldState<P,Q>& s) = 0;
    
    //! Leap with Q
    virtual void leapQ(const Real& dt, AbsFieldState<P,Q>& s) = 0;



    //! Get the trajectory length
    virtual const Real getTrajLength(void) const = 0;

    //! Get the step size 
    virtual const Real getStepSize(void) const = 0;
  };


}; // End namespace Chroma

#endif
