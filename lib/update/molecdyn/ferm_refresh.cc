#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "ferm_refresh.h"

using namespace QDP;

template<typename Hamiltonian>
void PseudoFermionHeatbath_t(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
			                      multi1d<LatticeColorMatrix>, 
			                      LatticeFermion >& mc_state,
			                      const Hamiltonian& H) 
{

  // Get at the gauge field
  const multi1d<LatticeColorMatrix>& u = mc_state.getQ();
  
  // Heatbath all the fields

  for(int pf = 0; pf < mc_state.getPhi().size(); pf++) { 
    
    // Get at the ferion action for piece i
    const FermionAction<LatticeFermion>& S_f = H.getFermAct(pf);
    
    // Create a Connect State, apply fermionic boundaries
    Handle< const ConnectState > f_state(S_f.createState(mc_state.getQ()));
    
    // Create a linear operator
    Handle< const LinearOperator<LatticeFermion> > D(S_f.linOp(f_state));
    
    LatticeFermion eta=zero;
    
    // Fill the eta field with gaussian noise
    gaussian(eta, D->subset());
    
    // Temporary: Move to correct normalisation
    eta *= sqrt(0.5);
    
    // Now HIT IT with the ROCK!!!! (Or in this case M^{dagger})
    (*D)((mc_state.getPhi())[pf], eta, MINUS);
  }				    
}

// For Inexact (MD) Hamiltonians (eg in R algorithm)
void PseudoFermionHeatbath(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
                                            multi1d<LatticeColorMatrix>, 
                                            LatticeFermion >& mc_state,
			   const AbsFermHamiltonian< 
                                            multi1d<LatticeColorMatrix>, 
                                            multi1d<LatticeColorMatrix>, 
                                            LatticeFermion >& H) 
{ 			  
  PseudoFermionHeatbath_t(mc_state,H);
}
