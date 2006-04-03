/* $Id: wldsduf_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/* This routine is specific to Wilson fermions! */

/* WlDsDuf -- computes the derivative of the fermionic action respect the  */
/*            link field */
/* u -- gauge field ( Read ) */

/*         |  dS      dS_f */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU */

/* psi -- [1./(M_dag*M)]*chi_  ( read ) */

void WlDsDuf(multi1d<LatticeColorMatrix>& ds_u,
	     const multi1d<LatticeColorMatrix>& u,
	     const LatticeFermion& psi)
{
  LINEAR_OPERATOR(A);

  LatticeColorMatrix utmp_1;
  LatticeFermion phi;
  LatticeFermion rho;
  LatticeFermion sigma;
  LatticeFermion ftmp_1;
  LatticeFermion ftmp_2;
  Double ddummy;
  Real dummy;
  int mu;
  int nu;
  int cb;
  
  START_CODE("subroutine");;
  
  
  /* Do the usual Wilson fermion dS_f/dU */
  CONSTRUCT_LINEAR_OPERATOR(A, lmpsi, u, KappaMD);

  /*  phi = M(u)*psi  */
  A (A, psi, phi, 1, PLUS);
    
  /* rho = Dslash(0<-1) * psi */
  dslash (u, psi, rho, PLUS, 1);
    
  /* phi = (KappaMD^2)*phi = -(KappaMD^2)*M*psi */
  dummy = -(KappaMD*KappaMD);
  phi = phi * dummy;
    
  /* sigma = Dslash_dag(0 <- 1) * phi */
  dslash (u, phi, sigma, MINUS, 1);
    
          
  for(mu = 0; mu < Nd; ++mu)
  {
    cb = 0;

    /* ftmp_2 = (gamma(mu))*psi */
    SPIN_PRODUCT(psi,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x )  */
    ftmp_2 -= psi;
    utmp_1 = -(shift(ftmp_2, cb, FORWARD, mu) * adj(sigma));

    /* ftmp_2 = (gamma(mu))*phi */
    SPIN_PRODUCT(phi,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x)  */
    ftmp_2 += phi;
    utmp_1 += shift(ftmp_2, cb, FORWARD, mu) * adj(rho);
    ds_u[mu][cb] += u[mu][cb] * utmp_1;
      
    cb = 1;

    /* ftmp_2 = (gamma(mu))*ftmp_1 */
    SPIN_PRODUCT(rho,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x )  */
    ftmp_2 -= rho;
    utmp_1 = -(shift(ftmp_2, cb, FORWARD, mu) * adj(phi));
      
    /* ftmp_2 = (gamma(mu))*sigma */
    SPIN_PRODUCT(sigma,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_1 = ftmp_1 + ftmp_2 = (1 + gamma(mu))*sigma( x + mu)  */
    ftmp_2 += sigma;
    utmp_1 += shift(ftmp_2, cb, FORWARD, mu) * adj(psi);
    ds_u[mu][cb] += u[mu][cb] * utmp_1;
      
  }

  FREE_LINEAR_OPERATOR(A);
  
  END_CODE("subroutine");;
}
