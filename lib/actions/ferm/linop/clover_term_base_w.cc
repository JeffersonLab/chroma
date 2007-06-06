// $Id: clover_term_base_w.cc,v 3.4 2007-06-06 22:10:19 bjoo Exp $
/*! \file
 *  \brief Clover term
 */

#include "actions/ferm/linop/clover_term_base_w.h"
#include "meas/glue/mesfield.h"

namespace Chroma 
{ 


  //! Return flops performed by the operator()
  unsigned long 
  CloverTermBase::nFlops() const {return 552;}


  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  void CloverTermBase::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			     const LatticeFermion& chi, const LatticeFermion& psi, 
			     enum PlusMinus isign) const
  {
    START_CODE();

    // base deriv resizes.
    // Even even checkerboard
    deriv(ds_u, chi, psi, isign,0);
    
    // Odd Odd checkerboard
    multi1d<LatticeColorMatrix> ds_tmp;
    deriv(ds_tmp, chi, psi, isign,1);
    
    ds_u += ds_tmp;
    
    END_CODE();
  }



  void CloverTermBase::deriv_loops(const int mu, const int nu, const int cb,
				   LatticeColorMatrix& ds_u,
				   const LatticeColorMatrix& Lambda,
				   const LatticeColorMatrix& Lambda_xplus_mu,
				   const LatticeColorMatrix& Lambda_xplus_nu,
				   const LatticeColorMatrix& Lambda_xplus_muplusnu) const
  {
    START_CODE();

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

    LatticeColorMatrix staple_for;
    LatticeColorMatrix staple_back;

    LatticeColorMatrix u_nu_for_mu = shift(u[nu],FORWARD, mu); // Can reuse these later
    LatticeColorMatrix u_mu_for_nu = shift(u[mu],FORWARD, nu);

    LatticeColorMatrix u_tmp1, u_tmp2, u_tmp3;
    LatticeColorMatrix u_tmp4,u_tmp5;
    //    LatticeColorMatrix Lambda_xplus_nu = shift(Lambda,FORWARD,nu);
    // LatticeColorMatrix Lambda_xplus_mu = shift(Lambda,FORWARD,mu);
    // LatticeColorMatrix Lambda_xplus_muplusnu = shift(Lambda_xplus_mu, FORWARD,nu);

    LatticeColorMatrix ds_tmp;

    //   u_tmp1 =   <-------
    //              |
    //              |                       ON ALL CHECKERBOARDS
    //              |                       (Because it's used in staples)  
    //              V
    u_tmp1 = adj(u_mu_for_nu)*adj(u[nu]);

    //  Staple_for =  <--------
    //                |       ^
    //                |       |             ON ALL CHECKERBOARDS
    //                |       |             
    //                V       |
    staple_for = u_nu_for_mu*u_tmp1;


    //
    //            ^
    //  u_tmp2 =  |                         ON ALL CHECKBERBOARDS
    //            |                         (Because it's used in the staples)
    //            |  
    //            <-------
    u_tmp2 = adj(u[mu])*u[nu];

    //  Staple_back =  ^       |
    //                 |       |            ON ALL CHECKERBOARDS
    //                 |       |
    //                 <------ V
    //                      
    staple_back = adj(u_nu_for_mu)*u_tmp2;


    // Now compute the terms of the force:
    // 
    // Altogether 8 terms. 4 Upwards with + sign, and 4 Downwards with - sign
    //                     4 terms use staples and 4 don't

    // NON STAPLE TERMS FIRST:

    // 1)
    //
    //    <-------  X (CB)                      <--------
    //    |         ^                           |
    //    |         |        re  use  u_tmp1 =  | 
    //    V         |                           V
    //   CB       1-CB       
    u_tmp3[rb[cb]] = u_nu_for_mu*Lambda_xplus_muplusnu;
    ds_u[rb[cb]] = u_tmp3*u_tmp1;
    
    // 2)
    //
    //  CB
    //    X <------ 
    //    |        ^       re use u[nu](x+mu) = u_nu_for_mu
    //    |        |       re use u[mu](x+nu) = u_mu_for_nu
    //    V        |
    //    1-CB    CB        
    u_tmp3[rb[1-cb]] = Lambda_xplus_nu*adj(u[nu]);
    u_tmp4[rb[1-cb]] = u_nu_for_mu*adj(u_mu_for_nu);
    ds_u[rb[1-cb]] = u_tmp4 * u_tmp3;

    // Terms 3) and 4)
    //
    // These last two can be done on the other checkerboard and then shifted together. at the very end...
    //
    //  CB      1-CB          
    //    ^       |   
    //    |       |     
    //    |       V
    //    <-------X CB

    // Compute on CB:     | 
    //                    |
    //         u_tmp1 =   |
    //                    V
    //                    X (CB)
    // u_tmp1[rb[cb]] = adj(u[nu])*Lambda;

    // 3)
    //
    //  Compunte with u_tmp2:     ^           |
    //                            |           | 
    //                  u_tmp2 =  |           | = u_tmp1
    //                            |           V
    //                     (1-CB) <--------   X CB
    u_tmp4[rb[1-cb]] = adj(u_nu_for_mu)*Lambda_xplus_mu;
    ds_tmp[rb[1-cb]] = u_tmp4*u_tmp2;
    
    // 4)
    //
    //  1-CB      CB
    //   ^        |
    //   |        |        reuse = u[nu](x+mu) = u_nu_for_mu
    //   |        |
    //   X <----- V 1-CB
    //   CB
    u_tmp4[rb[cb]] = Lambda*u[nu];
    u_tmp5[rb[cb]] = adj(u[mu])*u_tmp4;
    ds_tmp[rb[cb]] = adj(u_nu_for_mu)*u_tmp5;
    
    //  ds_tmp now holds the last 2 terms, one on each of its checkerboards, but Now I need
    //  to shift them both together onto ds_u
    //  I'll keep them in ds_tmp right, bearing in mind I'll need to bring
    //  them in with a -ve contribution...


    // STAPLE TERMS:
    
    // 5)
    //
    //    <------- 
    //    |        ^
    //    |        |     use computed staple
    //    V        |
    //    x        
    //   CB       1-CB  
    ds_u[rb[cb]] += staple_for*Lambda;
    
    // 6)
    //
    //    <------- 
    //    |        ^
    //    |        |    re use computed staple 
    //    V        |
    //   1-CB      X CB	  

    // lambda_tmp = shift(Lambda, FORWARD, mu);

    ds_u[rb[1-cb]] += Lambda_xplus_mu*staple_for;

    // 7)
    //
    //   CB      1-CB
    //  X         
    //    ^       |
    //    |       |   re use computed staple 
    //    |       |
    //    <-------V 
    //
    // This will need to be shifted back in nu -- add to ds_tmp
    // has -ve contribution, so makes a +ve contribution to ds_tmp
    ds_tmp[rb[1-cb]] += staple_back*Lambda_xplus_nu;


    // 8)
    //
    //  1-CB      X CB
    //    ^       |
    //    |       |  reuse computed staple 
    //    |       |
    //    <-------V
    //
    // This will need to be shifted back in nu -- add to ds_tmp
    // has a -ve contribution so makes a +ve contribution to ds_tmp
    ds_tmp[rb[cb]] += Lambda_xplus_muplusnu * staple_back;

    // Shift and accumulate ds_tmp with a -ve sign - both checkerboards
    u_tmp4 = shift(ds_tmp, BACKWARD, nu);
    ds_u -= u_tmp4;

    END_CODE();
  }


  //! Take deriv of D
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  void CloverTermBase::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			     const LatticeFermion& chi, const LatticeFermion& psi, 
			     enum PlusMinus isign, int cb) const
  {
    START_CODE();


    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    ds_u = zero;

    // Get the links
    const multi1d<LatticeColorMatrix>& u = getU();

    // spin trace gamma_mu gamma_nu X outer Y_dag for each mu, nu
    // (mu > nu)

    // Total number of nu > mu combinations
    int Ncomb = Nd*(Nd-1)/2;

    // spin trace gamma_mu gamma_nu X outer Y_dag for each mu, nu
    // (mu > nu)
    multi1d<LatticeColorMatrix> s_xy_dag(Ncomb);
    
    // shift of spinTraceOuterProd(x,y)  FORWARD mu
    multi1d<LatticeColorMatrix> s_xy_dag_pm(Ncomb);

    // shift of spinTraceOuterProd(x,y)  FORWARD nu
    multi1d<LatticeColorMatrix> s_xy_dag_pn(Ncomb);

    // shift of spinTraceOuterProd(x,y)  FORWARD nu + mu 
    multi1d<LatticeColorMatrix> s_xy_dag_pm_pn(Ncomb);

    // Make a compact index table
    multi2d<int> index(Nd,Nd);
    { 
      int foo = 0;
      for(int mu=0; mu < Nd; mu++) { 
	for(int nu=mu+1; nu < Nd; nu++) { 
	  index(mu,nu) = foo;  // Nu mu and mu nu are the same
	  index(nu,mu) = foo;  // up to a -ve sign from permuting
	  foo++;               // gamma_mu and gamma_nu
	}
      }
    }

    // Compute the actual spin traces
    for(int mu=0; mu < Nd; mu++) {
      for(int nu=mu+1; nu < Nd; nu++) { 

	Real factor = (Real(-1)/Real(8))*getCloverCoeff(mu,nu);
	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}
	LatticeFermion ferm_tmp = Gamma(mu_nu_index)*psi;
	int idx = index(mu,nu);
	s_xy_dag[idx] = 
	  factor*traceSpin( outerProduct(ferm_tmp,chi));
	s_xy_dag_pm[idx] = shift(s_xy_dag[idx], FORWARD, mu);
	s_xy_dag_pn[idx] = shift(s_xy_dag[idx], FORWARD, nu);
	s_xy_dag_pm_pn[idx] = shift(s_xy_dag_pm[idx], FORWARD, nu);
      }
    }

    // Now compute the insertions
    for(int mu=0; mu < Nd; mu++) 
    { 

      // I am only computing one checkerboard and intentionally
      // zeroing the other. This means I initialise both checkerboards
      // to zero since I will accumulate into the desired checkerboard
      // during the loops over mu and nu
      // ds_u[mu]  = zero;
      

      // Now we loop over nu and we build up 
      //
      // C1+, C2+, C3+, C4+ and C1- C2- C3- C4-
      //
      // These are sums over nu != mu but using symmetry we can 
      // write them as 2 sum nu > mu i sigma_mu F_munu
      //
      // 
      for(int nu = mu+1; nu < Nd; nu++) 
      {
	
	LatticeColorMatrix ds_tmp=zero;
	int idx = index(mu,nu);

	ds_tmp=zero;
	// computes +ve contrib to F[mu]
	// from mu,nu clover leaf
	deriv_loops(mu, nu, cb, ds_tmp, s_xy_dag(idx),
		    s_xy_dag_pm(idx),
		    s_xy_dag_pn(idx),
		    s_xy_dag_pm_pn(idx));

	ds_u[mu] += ds_tmp;
	ds_tmp=zero;

	// computes -ve contribs to F[nu] (because gamma_mu gamma_nu <-> gamma_nu gamma_mu
	// from mu,nu clover leaf
	deriv_loops(nu, mu, cb, ds_tmp, s_xy_dag(idx),
		    s_xy_dag_pn(idx),
		    s_xy_dag_pm(idx),
		    s_xy_dag_pm_pn(idx));

	ds_u[nu] -= ds_tmp;

      }
    }


    // Clear out the deriv on any fixed links
    getFermBC().zero(ds_u);
    END_CODE();
  }


  //! Take deriv of D using Trace Log
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$  
   */
  void CloverTermBase::derivTrLn(multi1d<LatticeColorMatrix>& ds_u, 
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();
    
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    // Force in each direction
    for(int mu=0; mu < Nd; mu++) 
    { 
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


	  LatticeColorMatrix sigma_XY_dag_xpm = shift(sigma_XY_dag, FORWARD, mu);
	  LatticeColorMatrix sigma_XY_dag_xpn = shift(sigma_XY_dag, FORWARD, nu);
	  LatticeColorMatrix sigma_XY_dag_xpmpn = shift(sigma_XY_dag_xpm, FORWARD, nu);

	  LatticeColorMatrix ds_tmp; 
	  deriv_loops(mu, nu, cb, ds_tmp, sigma_XY_dag,
		      sigma_XY_dag_xpm,sigma_XY_dag_xpn,sigma_XY_dag_xpmpn );

	  ds_u[mu] += ds_tmp;

	}  // End if mu != nu
  
      } // End loop over nu

    } // end of loop over mu
    

    // Not sure this is needed here, but will be sure
    getFermBC().zero(ds_u);
    
    END_CODE();
  }
      

} // End Namespace Chroma
