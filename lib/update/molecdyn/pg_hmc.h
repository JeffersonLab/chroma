#ifndef pg_hmc_h
#define pg_hmc_h

#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "update/molecdyn/abs_hyb_int.h"
#include "update/molecdyn/abs_hmc.h" 

class PureGaugeHMCTraj : 
public AbsNoFermionHMCTraj< multi1d<LatticeColorMatrix>, 
                            multi1d<LatticeColorMatrix>,
                            AbsFieldState, 
                            ExactAbsHamiltonian,
                            AbsHamiltonian, 
                            AbsLatColMatHybInt
                          >
{
 public:
  ~PureGaugeHMCTraj() {};

  // Constructor
  //
  // Takes H_MC -- the MC Hamiltonian
  //       MD_Integrator -- the MD Integrator
  PureGaugeHMCTraj(const ExactAbsHamiltonian<
                             multi1d<LatticeColorMatrix>, 
                             multi1d<LatticeColorMatrix> >& H_MC_, 
		   const AbsLatColMatHybInt< 
                             AbsFieldState<multi1d<LatticeColorMatrix>,
                                           multi1d<LatticeColorMatrix> >,
                             AbsHamiltonian< multi1d<LatticeColorMatrix>,
                                             multi1d<LatticeColorMatrix> >  
                          >&  MD_integrator_) : 
    H_MC(H_MC_.clone()), MD_integrator(MD_integrator_.clone()) {
    traj_done = 0;
  }

  // Satisfy interfaces
  const ExactAbsHamiltonian<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& getMCHamiltonian(void) const {
    return *H_MC;
  }

  const AbsLatColMatHybInt< AbsFieldState<multi1d<LatticeColorMatrix>,
                                          multi1d<LatticeColorMatrix> >,
                            AbsHamiltonian< multi1d<LatticeColorMatrix>,
                                            multi1d<LatticeColorMatrix> > >& 
  getMDIntegrator(void) const { 
    return *MD_integrator;
  }


  const int& getTrajNum(void) const { 
    return traj_done;
  }

  int& getTrajNum(void) {
    return traj_done;
  }

  // Inherit operator() and acceptReject() 

private:
  int traj_done;
  Handle< const ExactAbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > H_MC;
  Handle< const AbsLatColMatHybInt<AbsFieldState<multi1d<LatticeColorMatrix>,
                                                 multi1d<LatticeColorMatrix> >,
                                   AbsHamiltonian< multi1d<LatticeColorMatrix>,
                                                   multi1d<LatticeColorMatrix> > > > MD_integrator;
		   
};



#endif
