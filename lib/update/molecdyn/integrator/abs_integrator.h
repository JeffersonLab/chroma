// -*- C++ -*-
// $Id: abs_integrator.h,v 3.7 2009-04-17 02:05:37 bjoo Exp $

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

    //! Reset any chronological predictors for the integrator
    virtual void resetPredictors(void) const = 0;
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

    //! Reset Integrators in this level and sub integrators
    virtual void resetPredictors(void) const {
      resetPredictorsThisLevel();  // This level 
      getSubIntegrator().resetPredictors(); // Recurse down
    }

  protected:
    //! Refresh fields in just this level
    virtual void refreshFieldsThisLevel(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const = 0;

    virtual void resetPredictorsThisLevel(void) const = 0;
  };
  

  //! New MD integrator interface 
  /*! @ingroup integrator */
  template<typename P, typename Q>
  class AbsMDIntegrator {
  public:

    //! Virtual destructor
    virtual ~AbsMDIntegrator(void) {}

    //! Do the trajectory for length trajLength 
    virtual void operator()(AbsFieldState<P,Q>&s, const Real& trajLength) const {
      // Default behaviour - get the subintegrator and integrate
      // the trajectory of length trajLength
      AbsComponentIntegrator<P,Q>& theIntegrator = getIntegrator();

      // This is a toplevel thing - not embeddable so here is the 
      // place to reset any predictors
      QDPIO::cout << "MD: Resetting Chrono Predictors at start of trajectory" << endl;
      theIntegrator.resetPredictors();

      // This is recursive so no further resets in here.
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


} // End namespace Chroma

#endif
