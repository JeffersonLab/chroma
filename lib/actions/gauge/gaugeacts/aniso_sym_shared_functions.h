#ifndef __aniso_sym_shared_functions_h__
#define __aniso_sym_shared_functions_h__

#include "chromabase.h"

// \brief Shared functions for anisotropic Symanzik improved actions
namespace Chroma { 
  namespace AnisoSym { 

    // Compute the contributions of a plaguette and rectangle in the mu.nu
    // plane to the MD force.
    //
    // Notes: since the mu,nu and nu,mu plaquettes are identical, the 
    // plaquette contribution will only be calculated if mu < nu.
    //
    // If one wants to omit the contributions of the rectangle that is 
    // 2 link long in the time dimension, one should set noTemporal2Link 
    // to true.
    
    // Suggested call sequence:
    // ds_tmp = zero;
    // for(int mu=0; mu < Nd; mu++) { 
    //    for(int nu=0; nu < Nd; nu++) { 
    //      if( mu == nu ) continue;
    //
    //      deriv_part( mu, nu, t_dir,
    //                  c_0[mu][nu], c_1[mu][nu], 
    //                  noTemporal2Link, 
    //                  ds_tmp,
    //                  u )
    //
    //    }
    // }
    //
    // for(int mu=0; mu < Nd; mu++) { 
    //     ds_u[mu] = u[mu]*ds_tmp[mu];
    // }
    void deriv_part(const int mu, 
		    const int nu, 
		    const int t_dir,
		    const Real& c_plaq_munu, 
		    const Real& c_rect_munu, 
		    const bool noTemporal2Link,
		    multi1d<LatticeColorMatrix>& ds_tmp,
		    const multi1d<LatticeColorMatrix>& u) ;


    // Computes plaquettes of  \mu x \nu plaquette and a 2 \mu x \nu 
    // rectangle. Returns a LatticeReal in lgimp.
    //
    //
    // Suggested call sequence:
    // LatticeReal lgimp=zero;
    // for(int mu=0; mu < Nd; mu++) { 
    //    for(int nu=0; nu < Nd; nu++) { 
    //      if( mu == nu ) continue;
    //
    //      S_part( mu, nu, t_dir,
    //                  c_0[mu][nu], c_1[mu][nu], 
    //                  noTemporal2Link, 
    //                  lgimp
    //                  u )
    //
    //    }
    // }
    //
    // S = Real(-1)/Real(Nc)*sum(lgimp);
    void  S_part(int mu, int nu, int t_dir,
		 Real c_plaq_munu,
		 Real c_rect_munu, 
		 bool noTemporal2Link,
		 LatticeReal& lgimp,
		 const multi1d<LatticeColorMatrix>& u);
   
  }
}

#endif
