// -*- C++ -*-
// $Id: abs_integrator.h,v 2.0 2005-09-25 21:04:41 edwards Exp $

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
      int  num_steps = getNumSteps();

      Real tau0 = getTrajLength();

      write(xml_out, "n_steps", num_steps);
      write(xml_out, "dt", dt);
      write(xml_out, "dtby2", dtby2);
      write(xml_out, "tau0", tau0);

      // Start time = 0
      Real t = Real(0);
      Real p_t = Real(0);

      push(xml_out, "MDSteps");

      // First half step by leapP
      push(xml_out, "elem");
      write(xml_out, "t", p_t);
      leapP(dtby2, s);
      p_t += dtby2;
      pop(xml_out); // pop("elem");
  
      // First full step by leapQ
      push(xml_out, "elem");
      write(xml_out, "t", t);
      leapQ(dt, s);
      t += dt;
      pop(xml_out);
	
      for(int step=0; step < num_steps-1; step++) { 
	// Combined 2 half steps for leapP
	push(xml_out, "elem");
	write(xml_out, "t", p_t);
	leapP(dt,s);
	p_t += dt;
	pop(xml_out);

	// Full step for leapQ
	push(xml_out, "elem");
	write(xml_out, "t", t);
	leapQ(dt,s);
	t += dt;
	pop(xml_out);

      }
       
      // Last half step for P
      push(xml_out, "elem");
      write(xml_out, "t",p_t);
      leapP(dtby2, s);
      p_t += dtby2;
      pop(xml_out);

      // Write out both times at the end of traj (should be the same)
      push(xml_out, "EndOfTraj");
      write(xml_out, "t", t);
      write(xml_out, "t_mom", p_t);
      pop(xml_out); // EndOfTraj;

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

    //! Get the number of timesteps
    virtual const int getNumSteps(void) const = 0;
  };


}; // End namespace Chroma

#endif
