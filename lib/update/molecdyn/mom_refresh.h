#ifndef MOM_REFRESH_H
#define MOM_REFRESH_H

#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"

using namespace QDP;

// Generic Templated version
template<typename FS, typename Hamiltonian>
void MomRefreshGaussian_t(FS& state, const Hamiltonian& H);

// For pure gauge fields states, inexact Hamiltonians
void MomRefreshGaussian(AbsFieldState< multi1d<LatticeColorMatrix>, 
			               multi1d<LatticeColorMatrix> >& state, 
			const AbsHamiltonian< 
                                       multi1d<LatticeColorMatrix>, 
			               multi1d<LatticeColorMatrix> >& H);

// For pure gauge field states, Exact Hamiltonians

void MomRefreshGaussian(AbsFieldState< multi1d<LatticeColorMatrix>, 
                                       multi1d<LatticeColorMatrix> >& state, 
                                       const ExactAbsHamiltonian<
                                              multi1d<LatticeColorMatrix>, 
                                              multi1d<LatticeColorMatrix> >& H);

// For Pseudofermionic field states, inxact Hamiltonians
void MomRefreshGaussian(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
                                         LatticeFermion>& state, 
			const AbsFermHamiltonian< multi1d<LatticeColorMatrix>, 
                                                  multi1d<LatticeColorMatrix>, 
                                                  LatticeFermion>& H);


// For Pseudofermionic Field states, Exact Hamiltonians
void MomRefreshGaussian(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
                                         LatticeFermion>& state, 
                        const ExactAbsFermHamiltonian< 
                                         multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
                                         LatticeFermion>& H);

#endif
