#ifndef hmc_classes_h
#define hmc_classes_h

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

template<class EHS, class HI, class FS>
class HMCTraj {
public: 
  
  // Virtual destructor
  virtual ~HMCTraj() {}
  
  // Get at the Exact Hamiltonian -- It would be great if I could
  // enforce this to be exact. Perhaps I should request that it be
  // an AbsExactHamiltonian<P,Q> 
  virtual const EHS& getMCHamiltonian(void) const = 0;
	   
  // Get at the HybridIntegrator -- Again I could ask for this
  // to come from a base class but this one may be a complex beastie
  // with monitoring and such so it is bes to template this
  virtual const HI& getMDIntegrator(void)  const = 0;

  // Traj accessor
  virtual const int& getTrajNum(void) const = 0;

  // Traj mutator
  virtual int& getTrajNum(void) = 0;

  
  // Do the HMC trajectory
  virtual void operator()(FS& mc_state,
			  const bool doAccept,
			  XMLWriter& monitor)
  {
    const HI& integrator=getMDIntegrator();
    const EHS& H=getMCHamiltonian();

    // HMC Algorithm.
    // 1) Refresh momenta
    H.refreshP(mc_state);
    
    // 2) SaveState
    Handle< FS > old_state(mc_state.clone());
    
    
    // Info for user:
    push(monitor, "HMCTraj");
    write(monitor, "TrajNum", getTrajNum());

    // 4) Integrate MD trajectory
    integrator(mc_state, monitor);

    // 5) Measure energy of the new state

    // 3) Measure energy of the old state
    Double KE_old, PE_old;
    H.mesE(*old_state, KE_old, PE_old);

    Double KE, PE;
    H.mesE(mc_state, KE, PE);

    Double DeltaKE = KE - KE_old;
    Double DeltaPE = PE - PE_old;
    Double DeltaH  = DeltaKE + DeltaPE;

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
