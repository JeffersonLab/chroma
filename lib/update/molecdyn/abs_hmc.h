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
#include "update/molecdyn/mom_refresh.h"
#include "update/molecdyn/ferm_refresh.h"
#include "update/molecdyn/global_metropolis_accrej.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "update/molecdyn/abs_hyb_int.h"

using namespace QDP;
using namespace std;


template<typename P, typename Q, typename Phi, 
	 template <typename, typename, typename> class FS, 
	 template <typename, typename, typename> class HMC,
	 template <typename, typename, typename> class HMD,
	 template <typename, typename> class HI>
class AbsFermionHMCTraj {
public: 
  
  // Virtual destructor
  virtual ~AbsFermionHMCTraj() {}
  
  // Get at the Exact Hamiltonian -- It would be great if I could
  // enforce this to be exact. Perhaps I should request that it be
  // an AbsExactHamiltonian<P,Q> 
  virtual const HMC<P,Q,Phi>& getMCHamiltonian(void) const = 0;
	   
  // Get at the HybridIntegrator -- Again I could ask for this
  // to come from a base class but this one may be a complex beastie
  // with monitoring and such so it is bes to template this
  virtual const HI<FS<P,Q,Phi>, HMD<P,Q,Phi> >& getMDIntegrator(void) const = 0;
  // Traj accessor
  virtual const int& getTrajNum(void) const = 0;

  // Traj mutator
  virtual int& getTrajNum(void) = 0;

  
  // Do the HMC trajectory
  virtual void operator()(FS<P,Q,Phi>& mc_state,
			  const bool doAccept,
			  XMLWriter& monitor)
  {
    const HI< FS<P,Q,Phi>, HMD<P,Q,Phi> >& MD = getMDIntegrator();
    const HMC<P,Q,Phi>& H_MC = getMCHamiltonian();

    // HMC Algorithm.
    // 1) Refresh momenta & Pseudofermions
    //
    refreshP(mc_state);
    refreshPhi(mc_state);

    // 2) SaveState
    Handle< FS<P,Q,Phi> > old_state(mc_state.clone());
    
    // Info for user:
    push(monitor, "HMCTraj");
    write(monitor, "TrajNum", getTrajNum());

    // 3) Integrate MD trajectory
    MD(mc_state, monitor);

    // 4) Measure energy of the old state
    Double KE_old, GE_old, FE_old;
    H_MC.mesE(*old_state, KE_old, GE_old, FE_old);

    // 5) Measure the energy of the new state
    Double KE, GE, FE;
    H_MC.mesE(mc_state, KE, GE, FE);

    // Work out energy differences
    Double DeltaKE = KE - KE_old;
    Double DeltaGE = GE - GE_old;
    Double DeltaFE = FE - FE_old;
    Double DeltaH  = DeltaKE + DeltaGE + DeltaFE;

    // Write output
    write(monitor, "start_KE", KE_old);
    write(monitor, "end_KE", KE);
    write(monitor, "start_GE", GE_old);
    write(monitor, "end_GE", GE);
    write(monitor, "start_FE", FE_old);
    write(monitor, "end_FE", FE_old);

    write(monitor, "DeltaKE", DeltaKE);
    write(monitor, "DeltaGE", DeltaGE);
    write(monitor, "DeltaFE", DeltaFE);
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

	// I am going to refresh these so maybe these copies 
	//  for the momenta and the pseudofermions
	// are not really needed...
	//
	// restore momenta
	mc_state.getP() = old_state->getP();
	//
	// restore pseudofermions
	for(int pf = 0; pf < mc_state.getPhi().size(); pf++) { 
	  mc_state.getPhi()[pf] = old_state->getPhi()[pf];
	}

      }
    }

    pop(monitor);

    // Increase trajectory count
    getTrajNum()++;

  }


  virtual void refreshP(FS<P,Q,Phi>& state) {
    const HMC<P,Q,Phi>& H_MC = getMCHamiltonian();
    MomRefreshGaussian_t< FS<P,Q,Phi>, HMC<P,Q,Phi> >(state, H_MC);
  }

  virtual void refreshPhi(FS<P,Q,Phi>& state) const = 0;

  bool acceptReject(const Double& DeltaH, XMLWriter& monitor) { 
    return globalMetropolisAcceptReject(DeltaH, monitor);
  }

      
};

class AbsLatColMatLatFermHMCTraj : 
  public AbsFermionHMCTraj< multi1d<LatticeColorMatrix>, 
                            multi1d<LatticeColorMatrix>,
                            LatticeFermion, 
                            AbsPFFieldState, 
                            ExactAbsFermHamiltonian, 
                            AbsFermHamiltonian,
                            AbsLatColMatHybInt>
{
public: 
  // Virtual destructor
  virtual ~AbsLatColMatLatFermHMCTraj() {}

  // Refresh Phi with the MC Hamiltonian
  virtual void refreshPhi(AbsPFFieldState<multi1d<LatticeColorMatrix>, 
			  multi1d<LatticeColorMatrix>, LatticeFermion >& state) const {

    const ExactAbsFermHamiltonian< multi1d<LatticeColorMatrix>, 
                                   multi1d<LatticeColorMatrix>,
                                   LatticeFermion >& H_MC = getMCHamiltonian();

    PseudoFermionHeatbath(state, H_MC);

  }
};


// Pure Gauge (No Fermion case)
template<typename P, typename Q, 
	 template<typename,typename> class FS, 
	 template<typename,typename> class HMC, 
	 template<typename,typename> class HMD,
	 template<typename,typename> class HI>
class AbsNoFermionHMCTraj {
public: 
  
  // Virtual destructor
  virtual ~AbsNoFermionHMCTraj() {}
  
  // Get at the Exact Hamiltonian -- It would be great if I could
  // enforce this to be exact. Perhaps I should request that it be
  // an AbsExactHamiltonian<P,Q> 
  virtual const HMC<P,Q>& getMCHamiltonian(void) const = 0;
	   
  // Get at the HybridIntegrator -- Again I could ask for this
  // to come from a base class but this one may be a complex beastie
  // with monitoring and such so it is bes to template this
  virtual const HI< FS<P,Q>, HMD<P,Q> >& getMDIntegrator(void)  const = 0;

  // Traj accessor
  virtual const int& getTrajNum(void) const = 0;

  // Traj mutator
  virtual int& getTrajNum(void) = 0;

  
  // Do the HMC trajectory
  virtual void operator()(FS<P,Q>& mc_state,
			  const bool doAccept,
			  XMLWriter& monitor)
  {
    const HI< FS<P,Q>, HMD<P,Q> >& MD = getMDIntegrator();
    const HMC<P,Q>& H_MC = getMCHamiltonian();

    // HMC Algorithm.
    // 1) Refresh momenta
    //
    // Overloaded function like QuarkProp
    MomRefreshGaussian(mc_state, H_MC);

    // 2) SaveState
    Handle< FS<P,Q> > old_state(mc_state.clone());
    
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


  virtual void refreshP(FS<P,Q>& state) {
    const HMC<P,Q>& H_MC = getMCHamiltonian();
    MomRefreshGaussian_t< FS<P,Q>, HMC<P,Q> >(state, H_MC);
  }


  virtual bool acceptReject(const Double& DeltaH, XMLWriter& monitor) const { 
    return globalMetropolisAcceptReject(DeltaH, monitor);
  }  
      
};




#endif
