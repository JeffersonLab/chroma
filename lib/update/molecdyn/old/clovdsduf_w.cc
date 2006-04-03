/* $Id: clovdsduf_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/* This routine is specific to Wilson fermions! */

/* ClovDsDuf -- computes the derivative of the fermionic action respect the  */
/*              link field */
/* u -- gauge field ( Read ) */

/*         |  dS      dS_f */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU */

/* psi -- [1./(M_dag*M)]*chi_  ( read ) */

void ClovDsDuf(multi1d<LatticeColorMatrix>& ds_u,
	       const multi1d<LatticeColorMatrix>& u,
	       LatticeFermion psi)
{
  LINEAR_OPERATOR(A);

  LATTICE_TRIANG(clov);
  LATTICE_TRIANG(invclov);
  LATTICE_FIELD_STRENGTH(f);
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
  int mu;
  int nu;
  int cb;
  int munu;
  int it1;
  int it2;
  
  START_CODE("subroutine");;
  
  
  /* Do the clover fermion dS_f/dU */
    
  /*+ */
  /* The Kappas used in dslash and clovms */
  /* Kappa_ds = Kappa_dslash = Kappa */
  /* Kappa_cl = Kappa_clovms = ClovCoeff * Kappa / u0^3 */
  /*- */
  Kappa_ds = KappaMD;
  Kappa_cl = ClovCoeff * KappaMD / (u0*u0*u0);
  if (AnisoP == YES)
    QDP_error_exit("anisotropy not supported");
    
          
  /* Calculate F(mu,nu) */
  MesField (u, f);
    
  /* Make the 'clover terms', respectively their inverse */
  makclov (f, clov, Kappa_cl, Kappa_cl, 0);
  chlclovms (clov, invclov, NO, ddummy);
  makclov (f, clov, Kappa_cl, Kappa_cl, 1);
    
      
          
  CONSTRUCT_LINEAR_OPERATOR(A, lclovmpsi, u, clov, invclov, Kappa_ds);

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
    
                
  for(mu = 0; mu < Nd; ++mu)
  {
    /*+ */
    /* First do the terms not involving any dA(cb)/dU */
    /*- */
    
    cb = 0;

    /* ftmp_2 = (gamma(mu))*psi */
    SPIN_PRODUCT(psi,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2(x) = -(psi(x) - ftmp_2(x)) = -(1 - gamma(mu))*psi( x + mu)  */
    ftmp_2 -= psi;
    /* utmp_4 = -ftmp_2(x+mu)*sigma_dag(x,mu) */
    utmp_4 = -(shift(ftmp_2, cb, FORWARD, mu) * adj(sigma));
      
    /* ftmp_2 = (gamma(mu))*phi */
    SPIN_PRODUCT(phi,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = phi + ftmp_2 = (1 + gamma(mu))*phi( x )  */
    ftmp_2 += phi;
    /* utmp_4 = ftmp_2(x+mu) * rho_dag */
    utmp_4 += shift(ftmp_2, cb, FORWARD, mu) * adj(rho);
    ds_u[mu][cb] += u[mu][cb] * utmp_4;

    cb = 1;

    /* ftmp_2 = (gamma(mu))*rho */
    SPIN_PRODUCT(rho,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = -(rho - ftmp_2) = -(1 - gamma(mu))*rho( x )  */
    ftmp_2 -= rho;
    utmp_4 = -(shift(ftmp_2, cb, FORWARD, mu) * adj(phi));
    ds_u[mu][cb] += u[mu][cb] * utmp_4;

    /* ftmp_2 = (gamma(mu))*sigma */
    SPIN_PRODUCT(sigma,(INTEGER_LSHIFT_FUNCTION(1,mu)),ftmp_2);
    /* ftmp_2 = sigma + ftmp_2 = (1 + gamma(mu))*sigma( x )  */
    ftmp_2 += sigma;
    utmp_4 = shift(ftmp_2, cb, FORWARD, mu) * adj(psi);
    ds_u[mu][cb] += u[mu][cb] * utmp_4;

    /*+ */
    /* Next do the terms involving a dA(cb)/dU */
    /*- */
                                
    for(nu = mu+1; nu < Nd; ++nu)
    {
      if ( nu == mu ) continue;

      munu = INTEGER_LSHIFT_FUNCTION(1,mu) + INTEGER_LSHIFT_FUNCTION(1,nu);

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
      SPIN_PRODUCT(psi, munu, ftmp_1);
      utmp_3 = ftmp_1 * adj(phi);
      utmp_6 = utmp_3 * dummy;

      /* tau = - Kappa_c/8 * sigma_{mu,nu} * rho * (-1) */
      /* (the last to correct for the factor included in sigma!) */
      dummy = - Kapcl8;
      SPIN_PRODUCT(rho, munu, ftmp_1);
      utmp_3 = ftmp_1 * adj(sigma);
      utmp_7 = utmp_3 * dummy;
	
      /*+ */
      /* Finally do terms with dA(0)/dU from the "tr(log(A(0)))" term */
      /*- */
      dummy = TO_REAL(2) * Kapcl8;
      triacntr (invclov, munu, utmp_1);
      utmp_8 = utmp_1 * dummy;
		
      for(it1=0; it1<2; it1++)
      {
	cb = 0;	  
	
	utmp_4(rb[(1-cb)]) = shift(utmp_7, FORWARD, mu);
	utmp_0 = utmp_6 * u[mu][1-cb];
	utmp_0 += u[mu][1-cb] * utmp_4;

	utmp_1 = adj(u[mu][1-cb]) * utmp_6;
	utmp_1 += utmp_4 * adj(u[mu][1-cb]);
	utmp_2 = u[mu][1-cb] * shift(utmp_8, (1-cb), FORWARD, mu);

	/* backward */
	utmp_5(rb[(1-cb)]) = shift(u[nu][cb], FORWARD, mu);
	utmp_3 = utmp_0;
	utmp_3 += utmp_2;
	utmp_4 = utmp_3 * utmp_5;
	staple_back = shift(adj[u[nu][1-cb]], cb, BACKWARD, nu) * shift(utmp_4, cb, BACKWARD, nu);
	  
	utmp_4 = adj(utmp_5) * utmp_1;
	staple_for = shift(utmp_4, cb, BACKWARD, nu) * shift(u[nu][1-cb], cb, BACKWARD, nu);

	utmp_4 = u[mu][1-cb] * utmp_5;
	staple_bla_back = shift(adj[u[nu][1-cb]], cb, BACKWARD, nu) * shift(utmp_4, cb, BACKWARD, nu);
	  
	/* forward */
	utmp_5(rb[cb]) = shift(u[nu][1-cb], FORWARD, mu);
	utmp_4 = shift(utmp_0, cb, FORWARD, nu) * adj(utmp_5);
	staple_back -= u[nu][cb] * utmp_4;	
	  
	utmp_1 -= adj(utmp_2);
	utmp_3 = utmp_5 * shift(utmp_1, cb, FORWARD, nu);
	staple_for -= utmp_3 * adj(u[nu][cb]);	

	utmp_4 = shift(u[mu][1-cb], cb, FORWARD, nu) * adj(utmp_5);
	staple_bla_for = u[nu][cb] * utmp_4;	

	  
	/* put in the final links */
	staple_for -= adj(staple_bla_for) * utmp_8;
	staple_back += utmp_8 * staple_bla_back;
	  
	staple_bla_back -= staple_bla_for;
	staple_for += adj(staple_bla_back) * utmp_7;		
	staple_back += utmp_7 * staple_bla_back;

	utmp_3(rb[cb]) = shift(utmp_6, FORWARD, mu);
	staple_for += utmp_3 * adj(staple_bla_back);
	staple_back += staple_bla_back * utmp_3;	
	  
	if(it1==0)
	{
	  ds_u[mu][cb] += u[mu][cb] * staple_for;
	  ds_u[mu][cb] += staple_back * adj(u[mu][cb]);
	}
	else
	{
	  ds_u[mu][cb] -= u[mu][cb] * staple_for;
	  ds_u[mu][cb] -= staple_back * adj(u[mu][cb]);
	}
	
	cb = 1;	  
	
	utmp_4(rb[(1-cb)]) = shift(utmp_6, FORWARD, mu);
	utmp_0 = utmp_7 * u[mu][1-cb];
	utmp_0 += u[mu][1-cb] * utmp_4;

	utmp_1 = adj(u[mu][1-cb]) * utmp_7;
	utmp_1 += utmp_4 * adj(u[mu][1-cb]);
	utmp_2 = utmp_8 * u[mu][1-cb];


	/* backward */
	utmp_5(rb[(1-cb)]) = shift(u[nu][cb], FORWARD, mu);
	utmp_3 = utmp_0;
	utmp_3 += utmp_2;
	utmp_4 = utmp_3 * utmp_5;
	staple_back = shift(adj[u[nu][1-cb]], cb, BACKWARD, nu) * shift(utmp_4, cb, BACKWARD, nu);
	
	utmp_4 = adj(utmp_5) * utmp_1;
	staple_for = shift(utmp_4, cb, BACKWARD, nu) * shift(u[nu][1-cb], cb, BACKWARD, nu);

	utmp_4 = u[mu][1-cb] * utmp_5;
	staple_bla_back = shift(adj[u[nu][1-cb]], cb, BACKWARD, nu) * shift(utmp_4, cb, BACKWARD, nu);
	  
	/* forward */
	utmp_5(rb[cb]) = shift(u[nu][1-cb], FORWARD, mu);
	utmp_4 = shift(utmp_0, cb, FORWARD, nu) * adj(utmp_5);
	staple_back -= u[nu][cb] * utmp_4;	

	
	utmp_1 -= adj(utmp_2);
	utmp_3 = utmp_5 * shift(utmp_1, cb, FORWARD, nu);
	staple_for -= utmp_3 * adj(u[nu][cb]);	

	utmp_4 = shift(u[mu][1-cb], cb, FORWARD, nu) * adj(utmp_5);
	staple_bla_for = u[nu][cb] * utmp_4;	

	/* put in the final links */
	utmp_3(rb[cb]) = shift(utmp_8, FORWARD, mu);
	staple_for -= utmp_3 * adj(staple_bla_for);
	staple_back += staple_bla_back * utmp_3;

	staple_bla_back -= staple_bla_for;
	staple_for += adj(staple_bla_back) * utmp_6;		
	staple_back += utmp_6 * staple_bla_back;
	  
	utmp_3(rb[cb]) = shift(utmp_7, FORWARD, mu);
	staple_for += utmp_3 * adj(staple_bla_back);
	staple_back += staple_bla_back * utmp_3;	
	  
	if(it1==0)
	{
	  ds_u[mu][cb] += u[mu][cb] * staple_for;
	  ds_u[mu][cb] += staple_back * adj(u[mu][cb]);
	}
	else
	{
	  ds_u[mu][cb] -= u[mu][cb] * staple_for;
	  ds_u[mu][cb] -= staple_back * adj(u[mu][cb]);
	}

	it2=mu;
	mu=nu;
	nu=it2;
      }      
	
    }
  }

  FREE_LINEAR_OPERATOR(A);
          
  END_CODE("subroutine");
}
