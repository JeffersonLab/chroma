#ifndef abs_hmc_h
#define abs_hmc_h

#include "chromabase.h"

// EHS is the exact Hamiltonian System
// HI  is the Hybrid Integrator
// FS  is the Field State

// One thing I can't do here, is to enforce that 
// EHS is actually an exact Hamiltonian System
// HI is  actually a Hybrid Integrator
// and FS is a Field State
// 
// However, potentially the HMC would even work for a fermionic system
// also, as long as I define PE to include the fermions as well, 
// however, this may be frustrated by the fact that I never include
// the fermionic update.

#include "update/field_state.h"
#include "update/molecdyn/abs_mom_refresh.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "update/molecdyn/abs_hyb_int.h"

using namespace QDP;
using namespace std;

template<typename P, typename Q>
class HMCTraj {
public: 
  
  // Virtual destructor
  virtual ~HMCTraj() {}
  
  // Get at the Exact Hamiltonian -- It would be great if I could
  // enforce this to be exact. Perhaps I should request that it be
  // an AbsExactHamiltonian<P,Q> 
  virtual const ExactAbsHamiltonian<P,Q>& getMCHamiltonian(void) const = 0;
	   
  // Get at the HybridIntegrator -- Again I could ask for this
  // to come from a base class but this one may be a complex beastie
  // with monitoring and such so it is bes to template this
  virtual const AbsHybInt<P,Q>& getMDIntegrator(void)  const = 0;

  // Get the momentum refreshment 
  virtual const AbsMomGaussianHeatbath<P,Q>& getMomRefresh(void) const = 0;

  // Traj accessor
  virtual const int& getTrajNum(void) const = 0;

  // Traj mutator
  virtual int& getTrajNum(void) = 0;

  
  // Do the HMC trajectory
  virtual void operator()(AbsFieldState<P,Q>& mc_state,
			  const bool doAccept,
			  XMLWriter& monitor)
  {
    const AbsHybInt<P,Q>& MD = getMDIntegrator();
    const ExactAbsHamiltonian<P,Q>& H_MC = getMCHamiltonian();
    const AbsMomGaussianHeatbath<P,Q>& MomRefresh = getMomRefresh();

    // HMC Algorithm.
    // 1) Refresh momenta
    MomRefresh(mc_state, H_MC);

    // 2) SaveState
    Handle< AbsFieldState<P,Q> > old_state(mc_state.clone());
    
    // Info for user:
    push(monitor, "HMCTraj");
    write(monitor, "TrajNum", getTrajNum());

    // 3) Integrate MD trajectory
    MD(mc_state, monitor);

    // 4) Measure energy of the old state
    Double KE_old, PE_old;
    H_MC.mesE(*old_state, KE_old, PE_old);

    // 5) Measure the energy of the new state
    Double KE, PE;
    H_MC.mesE(mc_state, KE, PE);

    // Work out energy differences
    Double DeltaKE = KE - KE_old;
    Double DeltaPE = PE - PE_old;
    Double DeltaH  = DeltaKE + DeltaPE;

    // Write output
    write(monitor, "start_KE", KE_old);
    write(monitor, "end_KE", KE);
    write(monitor, "start_PE", PE_old);
    write(monitor, "end_PE", PE);
    write(monitor, "DeltaKE", DeltaKE);
    write(monitor, "DeltaPE", DeltaPE);
    write(monitor, "DeltaH", DeltaH);
    write(monitor, "doAccept", doAccept);

    // If we intend to do an accept reject step
    // (ie we are not warming up)
    if( doAccept ) { 
     
      // Measure Acceptance
      bool acceptTestResult = acceptReject(DeltaH, monitor);
      write(monitor, "acceptTestResult", acceptTestResult);

      // If rejected restore fields
      // If accepted no need to do anything
      if ( ! acceptTestResult ) { 
	mc_state.getQ() = old_state->getQ();
	mc_state.getP() = old_state->getP();
      }
    }

    pop(monitor);

    // Increase trajectory count
    getTrajNum()++;

  }

  virtual bool acceptReject(const Double& DeltaH, XMLWriter& monitor) const { 
    // If deltaH is negative then always accept
    bool ret_val;

    push(monitor, "AcceptRejectTest");
    if ( toBool( DeltaH <= Double(0)) ) {
      ret_val = true;
      write(monitor, "AccProb", Double(1));

    }
    else {
      Double AccProb = exp(-DeltaH);
      Double uni_dev;
      random(uni_dev);
     
      write(monitor, "AccProb", AccProb);
      write(monitor, "random", uni_dev);
      
      if( toBool( uni_dev <= AccProb ) ) { 

	ret_val = true;

      }
      else {    

	ret_val = false;
      }
    }

    write(monitor, "AcceptState", ret_val);
    pop(monitor);

    return ret_val;
  }  
      
};


#endif
