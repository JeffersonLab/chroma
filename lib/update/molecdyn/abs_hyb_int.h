#ifndef abs_hyb_int_h
#define abs_hyb_int_h

#include "chromabase.h"
#include "update/field_state.h"

using namespace QDP;
using namespace std;


//! Abstract MD Integrator -- a'la Robert's previous attempts
template<typename P, typename Q, typename FS, typename H>
class AbsHybInt {
public:
  //virtual destructor
  virtual ~AbsHybInt() {};

  // Do the integration: The XMLBufferWriter is for monitoring
  virtual void operator()(FS& s, 
			  XMLWriter& mon_traj) const = 0;

  // Get step size
  virtual Real getStepSize(void) const = 0;

  // Get traj length
  virtual Real getTrajLength(void) const = 0;

  // Clone so we can copy through the base class
  virtual AbsHybInt<P,Q,FS,H>* clone(void) const = 0;
};

template<typename FS, typename H>
class AbsLatColMatHybInt : AbsHybInt< multi1d<LatticeColorMatrix>,
                                multi1d<LatticeColorMatrix>,
                                FS, H > {
public:
  //virtual destructor
  virtual ~AbsLatColMatHybInt() {};

  // Do the integration: The XMLBufferWriter is for monitoring
  virtual void operator()(FS& s, 
			  XMLWriter& mon_traj) const = 0;


  // Get step size
  virtual Real getStepSize(void) const = 0;

  // Get traj length
  virtual Real getTrajLength(void) const = 0;

  // Clone so we can copy through the base class
  virtual AbsLatColMatHybInt<FS,H>* clone(void) const = 0;
};

template<typename FS, typename H>
class AbsLatColMatPQPLeapFrog : public AbsLatColMatHybInt< FS, H > {
public:
  //virtual destructor
  virtual ~AbsLatColMatPQPLeapFrog() {};

  virtual const AbsLatColMatSympUpdates<FS,H>& getSympUpdates(void) const = 0;

  // Do the integration: The XMLBufferWriter is for monitoring
  virtual void operator()(FS& s, 
			  XMLWriter& mon_traj) const {

    
    const AbsLatColMatSympUpdates<FS,H>& leaps = getSympUpdates();
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
      if( toBool( fabs(tau0 - t) <  dtby2  ) ) {
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
  
  virtual AbsLatColMatPQPLeapFrog<FS,H>* clone(void) const = 0;
};



#if 0
//! Fermions
template<typename P, typename Q, typename Phi>
class AbsPFHybInt {
public:
  //virtual destructor
  virtual ~AbsPFHybInt() {};


  // Do the integration: The XMLBufferWriter is for monitoring
  virtual void operator()(AbsPFFieldState<P,Q,Phi>& s, 
			  XMLWriter& mon_traj) const = 0;

  // Get step size
  virtual Real getStepSize(void) const = 0;

  // Get traj length
  //  virtual Real getTrajLength(void) const = 0;
  virtual Real getTrajLength(void) const = 0;
  

};
#endif 

#endif
