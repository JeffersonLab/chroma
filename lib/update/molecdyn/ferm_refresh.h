#ifndef FERM_REFRESH_H
#define FERM_REFRESH_H

#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"

using namespace QDP;

template<typename Hamiltonian>
void PseudoFermionHeatbath_t(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
			                      multi1d<LatticeColorMatrix>, 
			                      LatticeFermion >& mc_state,
			     const Hamiltonian& H); 

// For Inexact (MD) Hamiltonians (eg in R algorithm)
// Exact Hamiltonian is OK as it inherits from inexact
void PseudoFermionHeatbath(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
                                            multi1d<LatticeColorMatrix>,
                                            LatticeFermion >& mc_state,
			   const AbsFermHamiltonian< 
			                    multi1d<LatticeColorMatrix>, 
                                            multi1d<LatticeColorMatrix>, 
			                    LatticeFermion >& H);

#endif
