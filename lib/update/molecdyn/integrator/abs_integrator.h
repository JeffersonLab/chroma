// -*- C++ -*-
// $Id: abs_integrator.h,v 3.4 2006-11-20 19:15:02 bjoo Exp $

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
  // 
  class NoSubIntegratorsException {
  public : 
    NoSubIntegratorsException() {}
  };
  
  //! MD integrator that can be used as a component for other integrators
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class AbsComponentIntegrator { 
  public:
    //! Virtual destructor
    virtual ~AbsComponentIntegrator(void) {} 

    //! Do an integration of length n*delta tau in n steps.
    virtual void operator()(AbsFieldState<P,Q>& s, 
			    const Real& traj_length) const = 0;

    //! Refresh fields in this level of the integrator (for R like algorithms)
    virtual void refreshFields(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const = 0;
  };

  //! MD component integrator that has a sub integrator (recursive)
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class AbsRecursiveIntegrator : public AbsComponentIntegrator<P,Q> { 
  public: 
    //! Virtual destructor 
    virtual ~AbsRecursiveIntegrator(void) {}

    //! Do an integration of lenght n*delta tau in n steps.
    virtual  void operator()(AbsFieldState<P,Q>& s, 
			     const Real& traj_length) const = 0;
			    

    //! Return the next level down integrator
    virtual AbsComponentIntegrator<P,Q>& getSubIntegrator() const = 0;

    //! Refresh fields in this level of the integrator and sub integrators.
    virtual void refreshFields(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const {
      refreshFieldsThisLevel(s); // Do the fields in this level
      getSubIntegrator().refreshFields(s); // Recurse down 
    }

  protected:
    //! Refresh fields in just this level
    virtual void refreshFieldsThisLevel(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const = 0;
  };
  

  //! New MD integrator interface 
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class AbsMDIntegratorNew {
  public:
    //! Virtual destructor
    virtual ~AbsMDIntegratorNew(void) {}

    //! Do the trajectory for length trajLength 
    virtual void operator()(AbsFieldState<P,Q>&s, const Real& trajLength) const {
      // Default behaviour - get the subintegrator and integrate
      // the trajectory of length trajLength
      AbsComponentIntegrator<P,Q>& theIntegrator = getIntegrator();
      theIntegrator(s, trajLength);
    }
    
    //! Refresh fields in the sub integrators (for R-like algorithms)
    virtual void refreshFields(AbsFieldState<P,Q>&s ) const { 
      getIntegrator().refreshFields(s); // Recursively refresh fields
    }
    
    //! Get the trajectory length
    virtual Real getTrajLength(void) const = 0;

    //! Copy equivalent fields into MD monomals before integration
    /*! It is up to the toplevel integrator to keep track of which 
        fields it needs to copy internally so that this function doesn't
	need its details exposed */
    virtual void copyFields(void) const = 0;
  private:

    //! Get the toplevel sub integrator
    virtual AbsComponentIntegrator<P,Q>& getIntegrator() const = 0;
    
  };

  //! MD integrator interface
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class  AbsMDIntegrator {
  public:

    // Virtual destructor
    virtual ~AbsMDIntegrator(void) {}

    // Do a trajectory -- that is all HMD/HMC needs to know
    virtual void operator()(AbsFieldState<P,Q>& s, const Real& traj_length) = 0;

    //! Get at the MD Hamiltonian
    //  This needs to be exposed to initialise the PF fields.
    virtual AbsHamiltonian<P,Q>& getHamiltonian(void)  = 0;

    virtual const Real getTrajLength() const = 0;
  }; // End class MD integrator


  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class PQPLeapfrogIntegrator : public AbsMDIntegrator<P,Q> {
  public:
    //! Virtual destructor
    virtual ~PQPLeapfrogIntegrator(void) {} 


    //! Perform a trajectory 
    virtual void operator()(AbsFieldState<P,Q>& s, const Real& traj_length) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      // Self encapsulation rule
      push(xml_out, "PQPLeapfrogIntegrator");

      Real dt = getStepSize(); 
      Real dtby2 = dt / Real(2);
      int  num_steps = getNumSteps();

      // Real tau0 = getTrajLength();

      write(xml_out, "n_steps", num_steps);
      write(xml_out, "dt", dt);
      write(xml_out, "dtby2", dtby2);
      write(xml_out, "tau0", traj_length);

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
    
      END_CODE();
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
    virtual int getNumSteps(void) const = 0;
  };

} // End namespace Chroma

#endif
