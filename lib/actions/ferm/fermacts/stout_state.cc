// -*- C++ -*-
// $Id: stout_state.cc,v 2.1 2005-09-27 21:16:19 bjoo Exp $
/*! @file 
 *  @brief Connection State for Stout state (.cpp file)
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/stout_state.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"



namespace Chroma { 


  void StoutConnectState::smear_links(const multi1d<LatticeColorMatrix>& current, 
				      multi1d<LatticeColorMatrix>& next)
  {
    START_CODE();

    // Construct and add the staples, where smear_this_dirP[mu] is true
    LatticeColorMatrix u_staple = 0;
    for(int mu = 0; mu < Nd; mu++) { 

      for(int nu = 0; nu < Nd; ++nu) {
	
	if( nu != mu && smear_in_this_dirP[nu] == true ) {

	  // Forward staple
	  u_staple += rho(mu,nu)* current[nu] * shift(current[mu], FORWARD, nu) * adj(shift(current[nu], FORWARD, mu));

	  // Backward staple
	  // tmp_1(x) = u_dag(x,nu)*u(x,mu)*u(x+mu,nu)
	  LatticeColorMatrix tmp_1 = rho(mu, nu)*adj(current[nu]) * current[mu] * shift(current[nu], FORWARD, mu);
	  
	  // u_staple(x) += tmp_1_dag(x-nu)
	  //             += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu)
	  u_staple += shift(tmp_1, BACKWARD, nu);
	}
      }

      // The proto smeared link
      LatticeColorMatrix u_tmp = u_staple * adj(current[mu]);
      
      // Take the trace-less anti-hermitian projection of the staple
      taproj(u_tmp);
      
      // Exactly exponentiate the Lie Algebra matrix
      // Now u_tmp = exp(iQ)
      expmat(u_tmp,EXP_EXACT);

      next[mu]=u_tmp*current[mu];
    }

    END_CODE();
  }




  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
				       const multi2d<Real>& sm_fact_,
				       const int n_smear_, 
				       const multi1d<bool>& smear_in_this_dirP_)
  {
    START_CODE();
    create(u_, sm_fact_, n_smear_, smear_in_this_dirP_);
    END_CODE();
  }


  //! Explicitly specify smearing factor tensor
  StoutState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
				const multi2d<Real>& sm_fact_,
				const int n_smear_) 
  {
    START_CODE();
    // Specify all smearings but no mask. Assume 
    // smearing desired in all directions
    multi1d<bool> smear_in_this_dirP_aux(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      smear_in_this_dirP_aux[mu] = true;
    }
    
    create(u_,. sm_fact_, n_smear_, smear_in_this_dirP_aux);
    END_CODE();
  }

     //! Construct isotropic smearing in all 4 directions
  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_, 
				       const Real& sm_fact_, 
				       const int   n_smear_) 
  { 
    START_CODE();

    multi2d<Real> sm_fact_array(Nd, Nd);
    multi1d<bool> smear_in_this_dirP_aux(Nd);
    
    // For each (mu,nu) set sm_fact_array(mu,nu)=sm_fact
    // (Isotropy). Since mu != nu ever, we set those
    // to zero for safety
    for(int mu=0; mu < Nd; mu++) { 
      
      for(int nu=0; nu < Nd; nu++) { 
	
	if( mu==nu ) { 
	  sm_fact_array[mu][nu] = 0;
	}
	else { 
	  sm_fact_array[mu][nu] = sm_fact_;
	}
	
      }
      
      // smearing in all 4 directions
      smear_in_this_dirP_aux[mu]=true;
    }
    
    // call the create
    create(u_, sm_fact_array, n_smear_, smear_in_this_dirP_aux);
    END_CODE();
  }

  //! Construct isotopic smearing in 3 directions 
  StoutConnectState::StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
				       const Real& sm_fact_, 
				       const int   n_smear,
				       const int   j_decay) 
  {
    START_CODE();
    multi2d<Real> sm_fact_array(Nd, Nd);
    multi1d<bool> smear_in_this_dirP_aux(Nd);
    
    // For each (mu,nu) set sm_fact_array(mu,nu)=sm_fact
    // (Isotropy). Since mu != nu ever, we set those
    // to zero for safety
    for(int mu=0; mu < Nd; mu++) { 
      
      for(int nu=0; nu < Nd; nu++) { 
	if( mu==nu ) { 
	  sm_fact_array[mu][nu] = 0;
	}
	else { 
	  sm_fact_array[mu][nu] = sm_fact_;
	}
      }
      
      // Mask out the j_decay direction
      if( mu != j_decay ) { 
	smear_in_this_dirP_aux[mu]=true;
      }
      else { 
	smear_in_this_dirP_aux[mu]=false;
      }
    }
      
    // call the create
    create(u_, sm_fact_array, n_smear_, smear_in_this_dirP_aux);
    END_CODE();
  }

  // create function
  void StoutConnectState::create(const multi1d<LatticeColorMatrix>& u_,
	      const multi2d<Real>& sm_fact_,
	      const int n_smear_, 
	      const multi1d<bool>& smear_in_this_dirP_) 
  { 
    START_CODE();

    // Copy smearing factors
    rho.resize(Nd, Nd);
    rho = sm_fact_;
    
    // set n_smear
    n_smear = n_smear_;

    // Copy the direction maske
    smear_in_this_dirP.resize(Nd);
    smear_in_this_dirP = smear_in_this_dirP_;

    // Allocate smeared links
    smeared_links.resize(n_smear + 1);
    for(int i=0; i <= n_smear; i++) { 
      smeared_links[i].resize(Nd);
    }
    

    // Copy thin links into smeared_links[0]
    for(mu=0; mu < Nd; mu++) { 
      smeared_links[0][mu] = u_[mu];
    }

    // Iterate up the smearings
    for(int i=1; i <= n_smear; i++) { 
      smear_links(smeared_links[i-1], smeared_links[i]);
    }
    END_CODE();
  }




}; // End namespace Chroma
