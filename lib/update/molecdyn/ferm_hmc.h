#ifndef ferm_hmc_h
#define ferm_hmc_h

#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "update/molecdyn/abs_hyb_int.h"
#include "update/molecdyn/abs_hmc.h" 

class GaugeFermHMCTraj : public AbsLatColMatLatFermHMCTraj
{
 public:
  ~GaugeFermHMCTraj() {}

  // Constructor
  //
  // Takes H_MC -- the MC Hamiltonian
  //       MD_Integrator -- the MD Integrator
  GaugeFermHMCTraj(const ExactAbsFermHamiltonian<
                             multi1d<LatticeColorMatrix>, 
                             multi1d<LatticeColorMatrix>,
                             LatticeFermion >& H_MC_, 
		   const AbsLatColMatHybInt< 
                             AbsPFFieldState<multi1d<LatticeColorMatrix>,
                                             multi1d<LatticeColorMatrix>,
                                             LatticeFermion >,
                             AbsFermHamiltonian< multi1d<LatticeColorMatrix>,
                                                 multi1d<LatticeColorMatrix>,
                                                 LatticeFermion >  
                          >&  MD_integrator_) : 
    H_MC(H_MC_.clone()), MD_integrator(MD_integrator_.clone()) {
    traj_done = 0;
  }

  // Satisfy interfaces
  const ExactAbsFermHamiltonian<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix>,LatticeFermion >& getMCHamiltonian(void) const {
    return *H_MC;
  }

  const AbsLatColMatHybInt< AbsPFFieldState<multi1d<LatticeColorMatrix>,
                                            multi1d<LatticeColorMatrix>,
                                            LatticeFermion >,
                            AbsFermHamiltonian< multi1d<LatticeColorMatrix>,
                                                multi1d<LatticeColorMatrix>,
                                                LatticeFermion > >& 
  getMDIntegrator(void) const { 
    return *MD_integrator;
  }


  const int& getTrajNum(void) const { 
    return traj_done;
  }

  int& getTrajNum(void) {
    return traj_done;
  }

  // Inherit rest
private:
  int traj_done;
  Handle< const ExactAbsFermHamiltonian<multi1d<LatticeColorMatrix>, 
                                        multi1d<LatticeColorMatrix>, 
                                        LatticeFermion> > H_MC;

  Handle< const AbsLatColMatHybInt<AbsPFFieldState<multi1d<LatticeColorMatrix>,
                                                 multi1d<LatticeColorMatrix>,
                                                 LatticeFermion >,
                                   AbsFermHamiltonian< multi1d<LatticeColorMatrix>,
                                                   multi1d<LatticeColorMatrix>,
                                                   LatticeFermion > 
                                   > 
        > MD_integrator;
		   

};



#endif
