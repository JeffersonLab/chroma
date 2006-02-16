// $Id: clover_term_base_w.cc,v 2.5 2006-02-16 02:24:46 bjoo Exp $
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
    // Even even checkerboard
    deriv(ds_u, chi, psi, isign,0);
    
    // Odd Odd checkerboard
    multi1d<LatticeColorMatrix> ds_tmp;
    deriv(ds_tmp, chi, psi, isign,1);
    
    ds_u += ds_tmp;

  }



  void CloverTermBase::deriv_loops(const int mu, const int nu, const int cb,
		   LatticeColorMatrix& ds_u,
		   const LatticeColorMatrix& Lambda) const
  {

    const multi1d<LatticeColorMatrix>& u = getU();

    // New thingie - now assume Lambda lives only on sites with checkerboard 
    // CB
    //            Lambda
    //   0           X           0           x = cb, O = 1-cb
    //
    //
    // Lambda                 Lambda 
    //   X           0           X
    //
    //
    //            Lambda 
    //   0           X           0
    //
    // So I can only construct 4 out of the 8 staples on the sites
    // that have CB and the OTHER 4 of the 8 staples on sites with 
    // 1-cb
    //
    
    // Sites with CB first:
    //
    // This is first pass, most dumbest most inefficient way
    //
    //    <------- 
    //    |        ^
    //    |        |
    //    V        |
    //    x        
    //   CB       1-CB  
	  
    ds_u[rb[cb]] = shift(u[nu], FORWARD, mu)
      *adj(shift(u[mu],FORWARD,nu))
      *adj(u[nu])
      *Lambda;

    //    <------- X (CB)
    //    |        ^
    //    |        |
    //    V        |
    //   CB       1-CB       
	  
    ds_u[rb[cb]] += shift(u[nu],FORWARD, mu)
      *shift( shift(Lambda, FORWARD, mu), FORWARD, nu)
      *adj(shift(u[mu], FORWARD, nu))
      *adj(u[nu]);
    

    //
    //  CB      1-CB         
    //    ^       |
    //    |       |
    //    |       V
    //    <-------X CB
    
    ds_u[rb[cb]] -= adj( shift(shift(u[nu], FORWARD, mu),BACKWARD,nu) )
      * shift(shift(Lambda, FORWARD, mu), BACKWARD, nu)
      * adj(shift(u[mu], BACKWARD, nu))
      * shift(u[nu], BACKWARD, nu);


    //
    //   CB      1-CB
    //  X         
    //    ^       |
    //    |       |
    //    |       |
    //    <-------V 
    
    ds_u[rb[cb]]  -= adj( shift(shift(u[nu], FORWARD, mu),BACKWARD,nu) )
      * adj(shift(u[mu], BACKWARD, nu))
      * shift(u[nu], BACKWARD, nu)
      * Lambda;


    // Now Sites with 1 - CB 
    //  CB
    //    X <------ 
    //    |        ^
    //    |        |
    //    V        |
    //    1-CB    CB        
	  
    ds_u[rb[1-cb]] = shift(u[nu], FORWARD, mu)
      *adj(shift(u[mu], FORWARD, nu))
      *shift(Lambda, FORWARD, nu)
      *adj(u[nu]);


    //    <------- 
    //    |        ^
    //    |        |
    //    V        |
    //   1-CB      X CB
	  
	  
    ds_u[rb[1-cb]] += shift(Lambda, FORWARD,mu)
      *shift(u[nu], FORWARD, mu)
      *adj(shift(u[mu],FORWARD,nu))
      *adj(u[nu]);

    //  1-CB      X CB
    //    ^       |
    //    |       |
    //    |       |
    //    <-------V
	  
    
    ds_u[rb[1-cb]] -= shift(Lambda, FORWARD, mu)
      * adj(shift( shift(u[nu], FORWARD, mu), BACKWARD, nu))
      * adj(shift( u[mu], BACKWARD, nu))
      * shift(u[nu], BACKWARD, nu);


    //  1-CB      CB
    //   ^        |
    //   |        | 
    //   |        |
    //   X <----- V 1-CB
    //   CB

    ds_u[rb[1-cb]] -= adj( shift(shift(u[nu], FORWARD, mu),BACKWARD,nu) )
      * adj(shift(u[mu], BACKWARD, nu))
      * shift(Lambda, BACKWARD, nu)
      * shift(u[nu], BACKWARD, nu);
    	
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

      // Get the links
      const multi1d<LatticeColorMatrix>& u = getU();

      // I am only computing one checkerboard and intentionally
      // zeroing the other. This means I initialise both checkerboards
      // to zero since I will accumulate into the desired checkerboard
      // during the loops over mu and nu
      ds_u[mu]  = zero;
      

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
	  Real factor = (Real(-1)/Real(8))*getCloverCoeff(mu,nu);
	  
	  // Account for gamma conventions
	  // This is because we can only represent gamma_i gamma_j with i < j
	  // if i > j we need to do it as - gamma_j gamma_i
	  if( nu < mu) { 
	    factor *= Real(-1);
	  }
	  
	  // Work out X Y^{\dagger} = Tr_spin gamma_mu gamma_nu F_mu_nu
	  LatticeFermion ferm_tmp = Gamma(mu_nu_index)*psi;
	  LatticeColorMatrix sigma_XY_dag = 
	    factor*traceSpin( outerProduct(ferm_tmp,chi));

	  LatticeColorMatrix ds_tmp;
	  deriv_loops(mu, nu, cb, ds_tmp, sigma_XY_dag);
	  ds_u[mu] += ds_tmp;
	} // End if mu != nu
	

      } // End loop over nu
      
    }

  }


  //! Take deriv of D using Trace Log
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of chi vector                  (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void CloverTermBase::derivTrLn(multi1d<LatticeColorMatrix>& ds_u, 
				 enum PlusMinus isign, int cb) const
  {
    
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    // Force in each direction
    for(int mu=0; mu < Nd; mu++) { 


      // I am only computing this checkerboard, so initialise this
      ds_u[mu] = zero;

      // Now we loop over nu and we build up 
      //
      // C1+, C2+, C3+, C4+ and C1- C2- C3- C4-
      //
      // These are sums over nu != mu but using symmetry we can 
      // write them as 2 sum nu > mu i sigma_mu F_munu
      //
      // 
 
      // WARNING WARNING: CHANGE THIS BACK TO nu = 0 WHEN TESTING IS DONE
      for(int nu = 0; nu < Nd; nu++) { 

	if ( mu != nu ) {

	  // Index 
	  int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}
	  // QDPIO::cout << "mu = " << mu << " nu= " << nu << endl;
	  // QDPIO::cout << "mu_nu_index = " << mu_nu_index << endl;

	  // The actual coefficient factor
	  Real factor = Real(-1)*getCloverCoeff(mu,nu)/Real(8);
	  
	  // Account for gamma conventions
	  // This is because we can only represent gamma_i gamma_j with i < j
	  // if i > j we need to do it as - gamma_j gamma_i
	  if( nu < mu) { 
	    factor *= Real(-1);
	  }
	  	
	  LatticeColorMatrix sigma_XY_dag=zero;

	  // Need Sigma On Both Checkerboards
	  triacntr(sigma_XY_dag, mu_nu_index, cb);

	  sigma_XY_dag[rb[cb]] *= factor;

	  LatticeColorMatrix ds_tmp; 
	  deriv_loops(mu, nu, cb, ds_tmp, sigma_XY_dag);
	  ds_u[mu] += ds_tmp;

	}  // End if mu != nu
  
      } // End loop over nu

    } // end of loop over mu
    
  }
      
    

} // End Namespace Chroma
