#ifndef abs_hyb_int_h
#define abs_hyb_int_h

#include "chromabase.h"
#include "update/field_state.h"

using namespace QDP;
using namespace std;


//! Abstract MD Integrator -- a'la Robert's previous attempts
template<typename P, typename Q>
class AbsHybInt {
public:
  //virtual destructor
  virtual ~AbsHybInt() {};

  // Do the integration: The XMLBufferWriter is for monitoring
  virtual void operator()(AbsFieldState<P,Q>& s, 
			  XMLWriter& mon_traj) const = 0;

  // Get step size
  virtual Real getStepSize(void) const = 0;

  // Get traj length
  virtual Real getTrajLength(void) const = 0;

  // Clone so we can copy through the base class
  virtual AbsHybInt* clone(void) const = 0;

};

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
