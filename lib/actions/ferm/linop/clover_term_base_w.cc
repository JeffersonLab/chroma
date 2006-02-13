// $Id: clover_term_base_w.cc,v 2.3 2006-02-13 01:20:34 bjoo Exp $
/*! \file
 *  \brief Clover term
 */

#include "actions/ferm/linop/clover_term_base_w.h"
#include "meas/glue/mesfield.h"

namespace Chroma 
{ 

  //! Return flops performed by the operator()
  unsigned long 
  CloverTermBase::nFlops() const {return 0;}     // NOTE: NEED TO FIGURE THIS OUT!!


  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void CloverTermBase::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			     const LatticeFermion& chi, const LatticeFermion& psi, 
			     enum PlusMinus isign) const
  {

    // base deriv resizes.
    deriv(ds_u, chi, psi, isign,0);
    
    // This is prudence -- even though I know we shouldn't resize
    // again
    multi1d<LatticeColorMatrix> ds_tmp;
    deriv(ds_tmp, chi, psi, isign,1);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu][rb[1]] = ds_tmp[mu];
    }
  }

  //! Take deriv of D
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of chi vector                  (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void CloverTermBase::deriv(multi1d<LatticeColorMatrix>& ds_u, 
				const LatticeFermion& chi, const LatticeFermion& psi, 
				enum PlusMinus isign, int cb) const
  {
    
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
   
    }

    // Force in each direction
    for(int mu=0; mu < Nd; mu++) { 


      // The components of the force corresponding to the upper
      // and lower parts of the diagram in Zs thesise
      multi1d<LatticeColorMatrix> C_plus(4);
      multi1d<LatticeColorMatrix> C_minus(4);

      // Get the links
      const multi1d<LatticeColorMatrix>& u = getU();

      // Init the clover leaf forces to zero;
      for(int i=0; i < 4; i++) { 
	C_plus[i] = zero;
	C_minus[i]= zero;
      }
      
      // I am only computing this checkerboard, so initialise this
      ds_u[mu][rb[cb]] = zero;

      // Now we loop over nu and we build up 
      //
      // C1+, C2+, C3+, C4+ and C1- C2- C3- C4-
      //
      // These are sums over nu != mu but using symmetry we can 
      // write them as 2 sum nu > mu i sigma_mu F_munu
      //
      // 
      for(int nu = 0; nu < Nd; nu++) { 

	if ( mu != nu ) {

	  // Index 
	  int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}

	  // The actual coefficient factor
	  Real factor = Real(-1)/Real(8);
	  
	  // Account for gamma conventions
	  // This is because we can only represent gamma_i gamma_j with i < j
	  // if i > j we need to do it as - gamma_j gamma_i
	  if( nu < mu) { 
	    factor *= Real(-1);
	  }
	  
	  // Work out X Y^{\dagger} = Tr_spin gamma_mu gamma_nu F_mu_nu
	  LatticeFermion ferm_tmp = Gamma(mu_nu_index)*psi;
	  LatticeColorMatrix sigma_XY_dag = 
	    factor*getCloverCoeff(mu,nu)*traceSpin( outerProduct(ferm_tmp,chi));



	  // This is first pass, most dumbest most inefficient way
	  //  C_plus_1
	  //
	  //    <------- 
	  //    |        ^
	  //    |        |
	  //    V        |
	  //    x        O
	  
	  
	  C_plus[0] += shift(sigma_XY_dag, FORWARD,mu)
	    *shift(u[nu], FORWARD, mu)
	    *adj(shift(u[mu],FORWARD,nu))
	    *adj(u[nu]);
	  
	  
	  //  C_plus_2
	  //    <------- O
	  //    |        ^
	  //    |        |
	  //    V        |
	  //    x        
	  
	  C_plus[1] += shift(u[nu],FORWARD, mu)
	    *shift( shift(sigma_XY_dag, FORWARD, mu), FORWARD, nu)
	    *adj(shift(u[mu], FORWARD, nu))
	    *adj(u[nu]);
	  
	  
	  //  C_plus_3
	  //    O <------- 
	  //    |        ^
	  //    |        |
	  //    V        |
	  //    x        
	  
	  C_plus[2] += shift(u[nu], FORWARD, mu)
	    *adj(shift(u[mu], FORWARD, nu))
	    *shift(sigma_XY_dag, FORWARD, nu)
	    *adj(u[nu]);
	  
	  
	  // C_plus_4
	  //
	  //     <------- 
	  //    |        ^
	  //    |        |
	  //    V        |
	  //  x 0        
	  
	  C_plus[3] += shift(u[nu], FORWARD, mu)
	    * adj(shift(u[mu], FORWARD, nu))
	    * adj(u[nu])
	    * sigma_XY_dag;
	  
	  
	  // C_minus_1 
	  //
	  //  x         O
	  //    ^       |
	  //    |       |
	  //    |       |
	  //    <-------V
	  
	  
	  C_minus[0] += shift(sigma_XY_dag, FORWARD, mu)
	    * adj(shift( shift(u[nu], FORWARD, mu), BACKWARD, nu))
	    * adj(shift( u[mu], BACKWARD, nu))
	    * shift(u[nu], BACKWARD, nu);
	  
	  
	  // C_minus_2 
	  //
	  //  x         
	  //    ^       |
	  //    |       |
	  //    |       |
	  //    <-------V O
	  
	  C_minus[1] += adj( shift(shift(u[nu], FORWARD, mu),BACKWARD,nu) )
	    * shift(shift(sigma_XY_dag, FORWARD, mu), BACKWARD, nu)
	    * adj(shift(u[mu], BACKWARD, nu))
	    * shift(u[nu], BACKWARD, nu);
	  
	  // C_minus_3
	  //
	  //  x         
	  //    ^       |
	  //    |       |
	  //    |       |
	  //  0 <-------V 
	  
	  C_minus[2] += adj( shift(shift(u[nu], FORWARD, mu),BACKWARD,nu) )
	    * adj(shift(u[mu], BACKWARD, nu))
	    * shift(sigma_XY_dag, BACKWARD, nu)
	    * shift(u[nu], BACKWARD, nu);
	  
	  
	  // C_minus_4
	  //
	  //  x         
	  //  O ^       |
	  //    |       |
	  //    |       |
	  //    <-------V 
	  
	  C_minus[3] += adj( shift(shift(u[nu], FORWARD, mu),BACKWARD,nu) )
	    * adj(shift(u[mu], BACKWARD, nu))
	    * shift(u[nu], BACKWARD, nu)
	    * sigma_XY_dag;
	

	} // End if mu != nu
  
      } // End loop over nu
      
      
      for(int i=0; i < 4; i++) { 
	
	ds_u[mu][rb[cb]] += C_plus[i];
	ds_u[mu][rb[cb]] -= C_minus[i];
	
      }
      
    }

  }
      
    

} // End Namespace Chroma
