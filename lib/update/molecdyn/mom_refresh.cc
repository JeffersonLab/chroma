
#include "chromabase.h"

#include "mom_refresh.h"
#include "util/gauge/taproj.h"

using namespace QDP;

template<typename FS, typename Hamiltonian>
void MomRefreshGaussian_t(FS& state, 
			  const Hamiltonian& H) 
{

  // Loop over direcsions
  for(int mu = 0; mu < Nd; mu++) {

    // Pull the gaussian noise
    gaussian(state.getP()[mu]);
    state.getP()[mu] *= sqrt(0.5);  // Gaussian Normalisation

    // Make traceless and antihermitian
    taproj(state.getP()[mu]);

  }
  
  // Zero out momenta on Zero Boundary
  H.applyPBoundary(state.getP());
}

// Wrapper for pure gauge field state, inexact Hamiltonian
void MomRefreshGaussian(AbsFieldState< multi1d<LatticeColorMatrix>, 
                                       multi1d<LatticeColorMatrix> >& state, 
			const AbsHamiltonian< multi1d<LatticeColorMatrix>, 
			                      multi1d<LatticeColorMatrix> >& H)
{
  MomRefreshGaussian_t(state, H);
}


// Wrapper for pure gauge field state, exact Hamiltonian
void MomRefreshGaussian(AbsFieldState< multi1d<LatticeColorMatrix>, 
                                       multi1d<LatticeColorMatrix> >& state, 
			const ExactAbsHamiltonian< 
			               multi1d<LatticeColorMatrix>, 
                                       multi1d<LatticeColorMatrix> >& H) 
{
  MomRefreshGaussian_t(state, H);
}


// Wrapper for pseudofermion field state, inexact Hamiltonian
void MomRefreshGaussian(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
                                         LatticeFermion>& state, 
			const AbsFermHamiltonian< 
			                 multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
                                         LatticeFermion>& H) 
{ 
  MomRefreshGaussian_t(state, H);
}


// Wrapper for pseudofermion  field state, exact Hamiltonian
void MomRefreshGaussian(AbsPFFieldState< multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
                                         LatticeFermion>& state, 
			const ExactAbsFermHamiltonian< 
			                 multi1d<LatticeColorMatrix>, 
                                         multi1d<LatticeColorMatrix>, 
			                 LatticeFermion>& H) 
{ 
  MomRefreshGaussian_t(state, H);
}

