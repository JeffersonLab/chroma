// $Id: clover_term_base_w.cc,v 2.2 2006-01-01 05:12:30 edwards Exp $
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


#if 1
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
    QDPIO::cerr << "Clover deriv: not implemented yet" << endl;
    QDP_abort(1);
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
    QDPIO::cerr << "Clover deriv: not implemented yet" << endl;
    QDP_abort(1);
  }

#else


  //! Derivative of even-odd preconditioned Clover dM/dU
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void 
  CloverTermBase::deriv(multi1d<LatticeColorMatrix>& ds_u,
			const LatticeFermion& chi, const LatticeFermion& psi, 
			enum PlusMinus isign) const
  {
    START_CODE();
  
    LatticeColorMatrix utmp_0;
    LatticeColorMatrix utmp_1;
    LatticeColorMatrix utmp_2;
    LatticeColorMatrix utmp_3;
    LatticeColorMatrix utmp_4;
    LatticeColorMatrix utmp_5;
    LatticeColorMatrix utmp_6;
    LatticeColorMatrix utmp_7;
    LatticeColorMatrix utmp_8;
    LatticeColorMatrix staple_for;
    LatticeColorMatrix staple_back;
    LatticeColorMatrix staple_bla_for;
    LatticeColorMatrix staple_bla_back;
    LatticeFermion phi;
    LatticeFermion rho;
    LatticeFermion sigma;
    LatticeFermion tau;
    LatticeFermion ftmp_1;
    LatticeFermion ftmp_2;
    Double ddummy;
    Real dummy;
    Real Kappa_ds;
    Real Kappa_cl;
    Real Kapcl8;
  
    /* Do the clover fermion dS_f/dU */

    /*  phi = M(u)*psi  */
    A (A, psi, phi, 1, PLUS);
    
    /* rho = A^-1(0) * Dslash(0<-1) * psi */
    dslash (u, psi, ftmp_1, PLUS, 1);
    ldumul (invclov, ftmp_1, rho);
    
    /* phi = (Kappa_ds^2)*phi = -(Kappa_ds^2)*M*psi */
    dummy = -(Kappa_ds*Kappa_ds);
    phi = phi * dummy;
    
    /* sigma = A^-1(0) * Dslash_dag(0 <- 1) * phi */
    /*       = -(Kappa_ds^2) * A^-1(0) * Dslash_dag(0 <- 1) * M * psi */
    dslash (u, phi, ftmp_1, MINUS, 1);
    ldumul (invclov, ftmp_1, sigma);
    
                
    for(int mu = 0; mu < Nd; ++mu)
    {
      /*+ */
      /* First do the terms not involving any dA(cb)/dU */
      /*- */
    
      cb = 0;

      /* ftmp_2 = (gamma(mu))*psi */
      ftmp_2 = Gamma(1 << mu) * psi;

      /* ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x + mu)  */
      ftmp_2 -= psi;
      /* utmp_4 = -ftmp_2(x+mu)*sigma_dag(x,mu) */
      utmp_4[rb[cb]] = -(shift(ftmp_2, FORWARD, mu) * adj(sigma));
      
      /* ftmp_2 = (gamma(mu))*phi */
      ftmp_2 = Gamma(1 << mu) * phi;

      /* ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x )  */
      ftmp_2 += phi;
      /* utmp_4 = ftmp_2(x+mu) * rho_dag */
      utmp_4[rb[cb]] += shift(ftmp_2, FORWARD, mu) * adj(rho);
      ds_u[mu][cb] += u[mu][cb] * utmp_4;

      cb = 1;

      /* ftmp_2 = (gamma(mu))*rho */
      ftmp_2 = Gamma(1 << mu) * rho;

      /* ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x )  */
      ftmp_2 -= rho;
      utmp_4[rb[cb]] = -(shift(ftmp_2, FORWARD, mu) * adj(phi));
      ds_u[mu][cb] += u[mu] * utmp_4;

      /* ftmp_2 = (gamma(mu))*sigma */
      ftmp_2 = Gamma(1 << mu) * sigma;

      /* ftmp_2 = sigma + ftmp_2 = (1 + gamma(mu))*sigma( x )  */
      ftmp_2 += sigma;
      utmp_4[rb[cb]] = shift(ftmp_2, FORWARD, mu) * adj(psi);
      ds_u[mu][cb] += u[mu] * utmp_4;

      /*+ */
      /* Next do the terms involving a dA(cb)/dU */
      /*- */
      for(int nu = mu+1; nu < Nd; ++nu)
      {
	if ( nu == mu ) continue;

	munu = (1 << mu) + (1 << nu);

	if ( nu < mu )
	  Kapcl8 = - Kappa_cl / TO_REAL(8);
	else
	  Kapcl8 = Kappa_cl / TO_REAL(8);

	/*+ */
	/* First do terms with dA(1)/dU */
	/*- */
	
	/* tau = - Kappa_c/8 * sigma_{mu,nu} * psi / (-Kappa_ds**2) */
	/* (the last to correct for the factor included in phi!) */
	dummy = - Kapcl8 / (Kappa_ds*Kappa_ds);
	ftmp_1 = Gamma(munu) * psi;
	utmp_3 = ftmp_1 * adj(phi);
	utmp_6 = utmp_3 * dummy;

	/* tau = - Kappa_c/8 * sigma_{mu,nu} * rho * (-1) */
	/* (the last to correct for the factor included in sigma!) */
	dummy = - Kapcl8;
	ftmp_1 = Gamma(munu) * rho;
	utmp_3 = ftmp_1 * adj(sigma);
	utmp_7 = utmp_3 * dummy;
	
	/*+ */
	/* Finally do terms with dA(0)/dU from the "tr(log(A(0)))" term */
	/*- */
	dummy = TO_REAL(2) * Kapcl8;
	triacntr (invclov, munu, utmp_1);
	utmp_8 = utmp_1 * dummy;
		
	for(int it1=0; it1<2; it1++)
	{
	  cb = 0;	  
	
	  utmp_4[rb[1-cb]] = shift(utmp_7, FORWARD, mu);
	  utmp_0[rb[1-cb]] = utmp_6 * u[mu];
	  utmp_0[rb[1-cb]] += u[mu] * utmp_4;

	  utmp_1[rb[1-cb]] = adj(u[mu]) * utmp_6;
	  utmp_1[rb[1-cb]] += utmp_4 * adj(u[mu]);
	  utmp_2[rb[1-cb]] = u[mu] * shift(utmp_8, FORWARD, mu);

	  /* backward */
	  utmp_5[rb[1-cb]] = shift(u[nu], FORWARD, mu);
	  utmp_3 = utmp_0;
	  utmp_3 += utmp_2;
	  utmp_4 = utmp_3 * utmp_5;
	  staple_back[rb[cb]] = shift(adj(u[nu])*utmp_4, BACKWARD, nu);
	  
	  utmp_4 = adj(utmp_5) * utmp_1;
	  staple_for[rb[cb]] = shift(utmp_4 * u[nu], BACKWARD, nu);

	  utmp_4[rb[1-cb]] = u[mu] * utmp_5;
	  staple_bla_back[rb[cb]] = shift(adj(u[nu])*utmp_4, BACKWARD, nu);
	  
	  /* forward */
	  utmp_5[rb[cb]] = shift(u[nu], FORWARD, mu);
	  utmp_4[rb[cb]] = shift(utmp_0, FORWARD, nu) * adj(utmp_5);
	  staple_back[rb[cb]] -= u[nu] * utmp_4;	
	  
	  utmp_1 -= adj(utmp_2);
	  utmp_3[rb[cb]] = utmp_5 * shift(utmp_1, FORWARD, nu);
	  staple_for[rb[cb]] -= utmp_3 * adj(u[nu]);	

	  utmp_4[rb[cb]] = shift(u[mu], FORWARD, nu) * adj(utmp_5);
	  staple_bla_for[rb[cb]] = u[nu] * utmp_4;	

	  
	  /* put in the final links */
	  staple_for -= adj(staple_bla_for) * utmp_8;
	  staple_back += utmp_8 * staple_bla_back;
	  
	  staple_bla_back -= staple_bla_for;
	  staple_for += adj(staple_bla_back) * utmp_7;		
	  staple_back += utmp_7 * staple_bla_back;

	  utmp_3[rb[cb]] = shift(utmp_6, FORWARD, mu);
	  staple_for += utmp_3 * adj(staple_bla_back);
	  staple_back += staple_bla_back * utmp_3;	
	  
	  if(it1==0)
	  {
	    ds_u[mu][rb[cb]] += u[mu] * staple_for;
	    ds_u[mu][rb[cb]] += staple_back * adj(u[mu]);
	  }
	  else
	  {
	    ds_u[mu][rb[cb]] -= u[mu] * staple_for;
	    ds_u[mu][rb[cb]] -= staple_back * adj(u[mu]);
	  }
	
	  cb = 1;	  
	
	  utmp_4[rb[1-cb]] = shift(utmp_6, FORWARD, mu);
	  utmp_0[rb[1-cb]] = utmp_7 * u[mu];
	  utmp_0[rb[1-cb]] += u[mu] * utmp_4;

	  utmp_1[rb[1-cb]] = adj(u[mu]) * utmp_7;
	  utmp_1[rb[1-cb]] += utmp_4 * adj(u[mu]);
	  utmp_2[rb[1-cb]] = utmp_8 * u[mu];


	  /* backward */
	  utmp_5[rb[1-cb]] = shift(u[nu], FORWARD, mu);
	  utmp_3 = utmp_0;
	  utmp_3 += utmp_2;
	  utmp_4 = utmp_3 * utmp_5;
	  staple_back[rb[cb]] = shift(adj(u[nu])*utmp_4, BACKWARD, nu);
	
	  utmp_4 = adj(utmp_5) * utmp_1;
	  staple_for[rb[cb]] = shift(utmp_4*u[nu], BACKWARD, nu);

	  utmp_4[rb[1-cb]] = u[mu] * utmp_5;
	  staple_bla_back[rb[cb]] = shift(adj(u[nu])*utmp_4, BACKWARD, nu);
	  
	  /* forward */
	  utmp_5[rb[cb]] = shift(u[nu], FORWARD, mu);
	  utmp_4[rb[cb]] = shift(utmp_0, FORWARD, nu) * adj(utmp_5);
	  staple_back[rb[cb]] -= u[nu] * utmp_4;	

	
	  utmp_1 -= adj(utmp_2);
	  utmp_3[rb[cb]] = utmp_5 * shift(utmp_1, FORWARD, nu);
	  staple_for[rb[cb]] -= utmp_3 * adj(u[nu]);	

	  utmp_4[rb[cb]] = shift(u[mu], FORWARD, nu) * adj(utmp_5);
	  staple_bla_for[rb[cb]] = u[nu] * utmp_4;	

	  /* put in the final links */
	  utmp_3[rb[cb]] = shift(utmp_8, FORWARD, mu);
	  staple_for -= utmp_3 * adj(staple_bla_for);
	  staple_back += staple_bla_back * utmp_3;

	  staple_bla_back -= staple_bla_for;
	  staple_for += adj(staple_bla_back) * utmp_6;		
	  staple_back += utmp_6 * staple_bla_back;
	  
	  utmp_3[rb[cb]] = shift(utmp_7, FORWARD, mu);
	  staple_for += utmp_3 * adj(staple_bla_back);
	  staple_back += staple_bla_back * utmp_3;	
	  
	  if(it1==0)
	  {
	    ds_u[mu][rb[cb]] += u[mu] * staple_for;
	    ds_u[mu][rb[cb]] += staple_back * adj(u[mu]);
	  }
	  else
	  {
	    ds_u[mu][rb[cb]] -= u[mu][cb] * staple_for;
	    ds_u[mu][rb[cb]] -= staple_back * adj(u[mu]);
	  }

	  it2=mu;
	  mu=nu;
	  nu=it2;
	}      
	
      }
    }

    END_CODE();
  }
#endif

} // End Namespace Chroma
