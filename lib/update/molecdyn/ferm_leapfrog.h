#ifndef fermion_leapfrog_h
#define fermion_leapfrog_h

#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_symp_updates.h"
#include "update/molecdyn/abs_hyb_int.h"

using namespace QDP;
using namespace std;

//! A Concrete Leapfrog for Pure Gauge Systems.
//  Templated on 

class GaugeFermPQPLeapFrog : public AbsLatColMatPQPLeapFrog<
  AbsPFFieldState< multi1d<LatticeColorMatrix>,
			   multi1d<LatticeColorMatrix>, LatticeFermion >,
  AbsFermHamiltonian< multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix>, LatticeFermion >
 >
{
 public: 
  // Virtual Destructor
  ~GaugeFermPQPLeapFrog() {}

  GaugeFermPQPLeapFrog(AbsGaugeFermSympUpdates& symp_updates_,
		       Real delta_tau_, 
		       Real tau_) : symp_updates(symp_updates_.clone()), delta_tau(delta_tau_), tau(tau_) {}
  
  // Get at the leap P and leap Q
  // through a Symplectic Updates step.
  // Use Base Class...
  const AbsGaugeFermSympUpdates& getSympUpdates(void) const {
    return *symp_updates;
  }

  // Get step size
  virtual Real getStepSize(void) const { return delta_tau; }

  // Get traj length
  //virtual Real getTrajLength(void) const { return tau; }
  virtual Real getTrajLength(void) const { return tau; }


  GaugeFermPQPLeapFrog(const GaugeFermPQPLeapFrog& f) : symp_updates(f.symp_updates->clone()), delta_tau(f.delta_tau), tau(f.tau) {}
  
  virtual GaugeFermPQPLeapFrog* clone(void) const {
    return new GaugeFermPQPLeapFrog(*this);
  }

protected:
  Handle<AbsGaugeFermSympUpdates> symp_updates;
  const Real delta_tau;
  const Real tau;
};
#endif
