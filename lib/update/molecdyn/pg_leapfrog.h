#ifndef pg_leapfrog_h
#define pg_leapfrog_h

#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_symp_updates.h"
#include "update/molecdyn/abs_hyb_int.h"

using namespace QDP;
using namespace std;

//! A Concrete Leapfrog for Pure Gauge Systems.
//  Templated on 

class PureGaugePQPLeapFrog : public AbsHybInt< multi1d<LatticeColorMatrix>, 
				multi1d<LatticeColorMatrix> >
{
 public: 
  // Virtual Destructor
  ~PureGaugePQPLeapFrog() {}

  PureGaugePQPLeapFrog(SymplecticUpdates<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& symp_updates_,
		       Real delta_tau_, 
		       Real tau_) : symp_updates(symp_updates_.clone()), delta_tau(delta_tau_), tau(tau_) {}
  
  // Get at the leap P and leap Q
  // through a Symplectic Updates step.
  // Use Base Class...
  const SymplecticUpdates<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& getSympUpdates(void) const {
    return *symp_updates;
  }

  // Get step size
  virtual Real getStepSize(void) const { return delta_tau; }

  // Get traj length
  //virtual Real getTrajLength(void) const { return tau; }
  virtual Real getTrajLength(void) const { return tau; }

  virtual void doReversedTraj(AbsFieldState<multi1d<LatticeColorMatrix>, 
			                    multi1d<LatticeColorMatrix> > &s,
			      XMLWriter& mon_traj) const
  {
    // DO Forward Traj
    (*this)(s, mon_traj);

    // Flip Momenta
    for(int mu=0; mu < Nd; mu++) { 
      s.getP()[mu] *= Double(-1);
    }

    // Do backward Traj
    (*this)(s, mon_traj);

    // Flip Momenta
    for(int mu=0; mu < Nd; mu++) { 
      s.getP()[mu] *= Double(-1);
    }
  }

  // This is the dumb one with no monitoring...
  // One can override...
  virtual void operator()(AbsFieldState<multi1d<LatticeColorMatrix>, 
			                multi1d<LatticeColorMatrix> > &s,
			  XMLWriter& mon_traj) const
  {
    const SymplecticUpdates<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& leaps = getSympUpdates();
    Real dt = getStepSize();
    Real dtby2 = dt / Real(2);

    Real tau0 = getTrajLength();
    Real t = Real(0);

    // Initialise endP to false
    bool endP =false;

    // Fist Leap P 
    leaps.leapP(s, dtby2);
    while(!endP) {
     
      // Then leap Q
      leaps.leapQ(s,dt);
      
      // Increment MD time
      t += dt;

      // Check if this was the last gauge update
      // Tricky because you cannot really check for t = tau
      // with any certainty due to potential FP equality testing issues
      // However if the time remaining is substantially less than dt
      // ie dtby2 is less than dt and is less than dt by any reasonable
      // epsilon then we finish
      if( toBool( tau  <  t + dt  ) ) {
	// Time left is less than dtby2
	// Finish with a half P leap and signal end
	leaps.leapP(s,dtby2);
	endP = true;

      }
      else { 
	// Time left is more than dtby2
	// Make a full P Leap
	leaps.leapP(s, dt);

      }
    }
  }


  PureGaugePQPLeapFrog(const PureGaugePQPLeapFrog& f) : symp_updates(f.symp_updates->clone()), delta_tau(f.delta_tau), tau(tau) {}
  
  virtual PureGaugePQPLeapFrog* clone(void) const {
    return new PureGaugePQPLeapFrog(*this);
  }

protected:
  Handle<SymplecticUpdates<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > symp_updates;
  Real delta_tau;
  Real tau;
};
#endif
