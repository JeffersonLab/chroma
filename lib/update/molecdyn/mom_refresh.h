#ifndef MOM_REFRESH_H
#define MOM_REFRESH_H

#include "chromabase.h"
#include "update/molecdyn/abs_mom_refresh.h"

// A Momentum update class. Works only for SU(3) matrices
class su3MomGaussianHeatbath :
  public AbsMomGaussianHeatbath< multi1d<LatticeColorMatrix>, 
                                   multi1d<LatticeColorMatrix> >
{
public: 
  virtual void operator()(AbsFieldState< multi1d<LatticeColorMatrix>, 
			                 multi1d<LatticeColorMatrix> >& state, 
			  const AbsHamiltonian< multi1d<LatticeColorMatrix>, 
			                        multi1d<LatticeColorMatrix> >& H) const {
    
    for(int mu = 0; mu < Nd; mu++) {
      gaussian(state.getP()[mu]);
      taproj(state.getP()[mu]);
      state.getP()[mu] *= sqrt(0.5);  // Gaussian Normalisation
    }
    
    // Zero out momenta on Zero Boundary
    H.applyPBoundary(state.getP());
  }
};

    

#endif
