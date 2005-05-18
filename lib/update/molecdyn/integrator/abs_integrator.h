// -*- C++ -*-
// $Id: abs_integrator.h,v 1.6 2005-05-18 18:30:12 edwards Exp $

/*! @file
 * @brief Integrators
 *
 * Intregators for HMC
 */

#ifndef ABS_INTEGRATOR_H
#define ABS_INTEGRATOR_H

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "io/xmllog_io.h"


namespace Chroma 
{

  //! MD integrator interface
  /*! @ingroup integrator */
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
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class PQPLeapfrogIntegrator : public AbsMDIntegrator<P,Q> {
  public:
    //! Virtual destructor
    virtual ~PQPLeapfrogIntegrator(void) {} 


    //! Perform a trajectory 
    virtual void operator()(AbsFieldState<P,Q>& s) {

      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      // Self encapsulation rule
      push(xml_out, "PQPLeapfrogIntegrator");

      Real dt = getStepSize(); 
      Real dtby2 = dt / Real(2);
      Real tau0 = getTrajLength();
      
      write(xml_out, "dt", dt);
      write(xml_out, "dtby2", dtby2);
      write(xml_out, "tau0", tau0);

      // Start time = 0
      Real t = Real(0);

      bool endP = false;

      push(xml_out, "MDSteps");

      // First half step by leapP
      push(xml_out, "elem");
      write(xml_out, "t", t);
      leapP(dtby2, s);
      pop(xml_out); // pop("elem");
      
      while(! endP ) { 
	push(xml_out, "elem");
	leapQ(dt, s);
	pop(xml_out);

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
	  push(xml_out, "elem");
	  write(xml_out, "t", t);
	  leapP(dtby2, s);
	  pop(xml_out);

	  endP = true;
	  
	}
	else {
	  push(xml_out, "elem");
	  write(xml_out, "t", t);
	  leapP(dt, s);
	  pop(xml_out);
	}
      } // end while
      pop(xml_out); // pop("MDSteps");
      pop(xml_out); // pop("PQPLeapfrogIntegrator")
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
