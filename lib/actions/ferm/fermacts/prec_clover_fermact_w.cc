// $Id: prec_clover_fermact_w.cc,v 1.5 2004-01-07 13:50:07 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_clover_linop_w.h"
#include "actions/ferm/fermacts/prec_clover_fermact_w.h"
#include "actions/ferm/linop/lmdagm.h"

//! Produce a linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecCloverFermAct::linOp(Handle<const ConnectState> state) const
{
  return new EvenOddPrecCloverLinOp(state->getLinks(),Mass,ClovCoeff,u0);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
EvenOddPrecCloverFermAct::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<LatticeFermion>(linOp(state));
}


//! Computes the derivative of the fermionic action respect to the link field
/*!
 *         |  dS      dS_f
 * ds_u -- | ----   + -----   ( Write )
 *         |  dU       dU
 *
 * psi -- [1./(M_dag*M)]*chi_  ( read ) 
 *
 * \param ds_u     result      ( Write )
 * \param state    gauge field ( Read )
 * \param psi      solution to linear system ( Read )
 */

void
EvenOddPrecCloverFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			       Handle<const ConnectState> state,
			       const LatticeFermion& psi) const
{
  START_CODE("EvenOddPrecCloverFermAct::dsdu");
  
  QDPIO::cerr << "EvenOddPrecCloverFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("EvenOddPrecCloverFermAct::dsdu");
}


#if 0

#error "NEEDS MORE CONVERSION"

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

/*# This routine is specific to Wilson fermions! */
/*# TRIACNTR - calculates */

/*#     Tr_D ( Gamma_mat L ) */

/*#  the trace over the Dirac indices for one of the 16 Gamma matrices */
/*#  and a hermitian color x spin matrix A, stored as a block diagonal */
/*#  complex lower triangular matrix L and a real diagonal diag_L. */

/*#  Here 0 <= mat <= 15 and */
/*#  if mat = mat_1 + mat_2 * 2 + mat_3 * 4 + mat_4 * 8 */

/*#  Gamma(mat) = gamma(1)^(mat_1) * gamma(2)^(mat_2) * gamma(3)^(mat_3) */
/*#             * gamma(4)^(mat_4) */

/*#  Further, in basis for the Gamma matrices used, A is of the form */

/*#      | A_0 |  0  | */
/*#  A = | --------- | */
/*#      |  0  | A_1 | */


/*# Arguments: */

/*#  L       -- lower triangular matrix		(Read) */
/*#  diag_L  -- diag(L)				(Read) */
/*#  mat     -- label of the Gamma matrix		(Read) */
/*#  B       -- the resulting SU(N) color matrix	(Write) */

void triacntr(LatticeColorMatrix& B,
	      const LATTICE_TRIANG& clov,
	      int mat)
{
  START_CODE("subroutine");;
  
  if ( mat < 0  ||  mat > 15 )
    QDP_error_exit("Gamma out of range", mat);
  
  switch( mat )
  {
  case 0:
    /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */
    /*#               0  1  0  0 */
    /*#               0  0  1  0 */
    /*#               0  0  0  1 */
    /*# From diagonal part */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp0;
	InnerReal lr_zero0;
	InnerReal lrtmp0;
	Complex sqm0;
	int i0;
	int j0;
	int elem_ij0;
	int elem_ijb0;
  
                          
	lr_zero0 = 0;
	FILL(sqm0, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i0 = 0; i0 < Nc; ++i0)
	    {
	      lrtmp0 = _DIAG_SP_clov[0][i0];
	      lrtmp0 += _DIAG_SP_clov[0][i0+Nc];
	      lrtmp0 += _DIAG_SP_clov[1][i0];
	      lrtmp0 += _DIAG_SP_clov[1][i0+Nc];
	      _SP_B[i0][i0] = cmplx(lrtmp0,lr_zero0);
	    }

	    /*# From lower triangular portion */
	    elem_ij0 = 0;
	    for(i0 = 1; i0 < Nc; ++i0)
	    {
	      elem_ijb0 = (i0+Nc)*(i0+Nc-1)/2 + Nc;
	      for(j0 = 0; j0 < i0; ++j0)
	      {
		lctmp0 = _OFFD_SP_clov[0][elem_ij0];
		lctmp0 += _OFFD_SP_clov[0][elem_ijb0];
		lctmp0 += _OFFD_SP_clov[1][elem_ij0];
		lctmp0 += _OFFD_SP_clov[1][elem_ijb0];

		_SP_B[j0][i0] = lctmp0;
		_SP_B[i0][j0] = adj(lctmp0);
	      
		elem_ij0 = elem_ij0 + 1;
		elem_ijb0 = elem_ijb0 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 3:
    /*# gamma(  12)  -i  0  0  0            # ( 0011 )  --> 3 */
    /*#               0  i  0  0 */
    /*#               0  0 -i  0 */
    /*#               0  0  0  i */
    /*# From diagonal part */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp3;
	InnerReal lr_zero3;
	InnerReal lrtmp3;
	Complex sqm3;
	int i3;
	int j3;
	int elem_ij3;
	int elem_ijb3;
  
                          
	lr_zero3 = 0;
	FILL(sqm3, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i3 = 0; i3 < Nc; ++i3)
	    {
	      lrtmp3 = _DIAG_SP_clov[0][i3+Nc];
	      lrtmp3 -= _DIAG_SP_clov[0][i3];
	      lrtmp3 -= _DIAG_SP_clov[1][i3];
	      lrtmp3 += _DIAG_SP_clov[1][i3+Nc];
	      _SP_B[i3][i3] = cmplx(lr_zero3,lrtmp3);
	    }
	
	    /*# From lower triangular portion */
	    elem_ij3 = 0;
	    for(i3 = 1; i3 < Nc; ++i3)
	    {
	      elem_ijb3 = (i3+Nc)*(i3+Nc-1)/2 + Nc;
	      for(j3 = 0; j3 < i3; ++j3)
	      {
		lctmp3 = _OFFD_SP_clov[0][elem_ijb3];
		lctmp3 -= _OFFD_SP_clov[0][elem_ij3];
		lctmp3 -= _OFFD_SP_clov[1][elem_ij3];
		lctmp3 += _OFFD_SP_clov[1][elem_ijb3];

		_SP_B[j3][i3] = lctmp3 * sqm3;
		_SP_B[i3][j3] = adj(lctmp3) * sqm3;
	    
		elem_ij3 = elem_ij3 + 1;
		elem_ijb3 = elem_ijb3 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 5:
    /*# gamma(  13)   0 -1  0  0            # ( 0101 )  --> 5 */
    /*#               1  0  0  0 */
    /*#               0  0  0 -1 */
    /*#               0  0  1  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp5;
	InnerReal lr_zero5;
	InnerReal lrtmp5;
	Complex sqm5;
	int i5;
	int j5;
	int elem_ij5;
	int elem_ji5;
  
                          
	lr_zero5 = 0;
	FILL(sqm5, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i5 = 0; i5 < Nc; ++i5)
	    {
	      elem_ij5 = (i5+Nc)*(i5+Nc-1)/2;
	      for(j5 = 0; j5 < Nc; ++j5)
	      {
		elem_ji5 = (j5+Nc)*(j5+Nc-1)/2 + i5;

		lctmp5 = adj(_OFFD_SP_clov[0][elem_ji5]);
		lctmp5 -= _OFFD_SP_clov[0][elem_ij5];
		lctmp5 += adj(_OFFD_SP_clov[1][elem_ji5]);
		lctmp5 -= _OFFD_SP_clov[1][elem_ij5];

		_SP_B[j5][i5] = lctmp5;

		elem_ij5 = elem_ij5 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 6:
    /*# gamma(  23)   0 -i  0  0            # ( 0110 )  --> 6 */
    /*#              -i  0  0  0 */
    /*#               0  0  0 -i */
    /*#               0  0 -i  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp6;
	InnerReal lr_zero6;
	InnerReal lrtmp6;
	Complex sqm6;
	int i6;
	int j6;
	int elem_ij6;
	int elem_ji6;
  
                          
	lr_zero6 = 0;
	FILL(sqm6, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i6 = 0; i6 < Nc; ++i6)
	    {
	      elem_ij6 = (i6+Nc)*(i6+Nc-1)/2;
	      for(j6 = 0; j6 < Nc; ++j6)
	      {
		elem_ji6 = (j6+Nc)*(j6+Nc-1)/2 + i6;

		lctmp6 = adj(_OFFD_SP_clov[0][elem_ji6]);
		lctmp6 += _OFFD_SP_clov[0][elem_ij6];
		lctmp6 += adj(_OFFD_SP_clov[1][elem_ji6]);
		lctmp6 += _OFFD_SP_clov[1][elem_ij6];

		_SP_B[j6][i6] = -(lctmp6 * sqm6);

		elem_ij6 = elem_ij6 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 9:
    /*# gamma(  14)   0  i  0  0            # ( 1001 )  --> 9 */
    /*#               i  0  0  0 */
    /*#               0  0  0 -i */
    /*#               0  0 -i  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp9;
	InnerReal lr_zero9;
	InnerReal lrtmp9;
	Complex sqm9;
	int i9;
	int j9;
	int elem_ij9;
	int elem_ji9;
  
                          
	lr_zero9 = 0;
	FILL(sqm9, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i9 = 0; i9 < Nc; ++i9)
	    {
	      elem_ij9 = (i9+Nc)*(i9+Nc-1)/2;
	      for(j9 = 0; j9 < Nc; ++j9)
	      {
		elem_ji9 = (j9+Nc)*(j9+Nc-1)/2 + i9;

		lctmp9 = adj(_OFFD_SP_clov[0][elem_ji9]);
		lctmp9 += _OFFD_SP_clov[0][elem_ij9];
		lctmp9 -= adj(_OFFD_SP_clov[1][elem_ji9]);
		lctmp9 -= _OFFD_SP_clov[1][elem_ij9];

		_SP_B[j9][i9] = lctmp9 * sqm9;

		elem_ij9 = elem_ij9 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 10:
    /*# gamma(  24)   0 -1  0  0            # ( 1010 )  --> 10 */
    /*#               1  0  0  0 */
    /*#               0  0  0  1 */
    /*#               0  0 -1  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp10;
	InnerReal lr_zero10;
	InnerReal lrtmp10;
	Complex sqm10;
	int i10;
	int j10;
	int elem_ij10;
	int elem_ji10;
  
                          
	lr_zero10 = 0;
	FILL(sqm10, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i10 = 0; i10 < Nc; ++i10)
	    {
	      elem_ij10 = (i10+Nc)*(i10+Nc-1)/2;
	      for(j10 = 0; j10 < Nc; ++j10)
	      {
		elem_ji10 = (j10+Nc)*(j10+Nc-1)/2 + i10;

		lctmp10 = adj(_OFFD_SP_clov[0][elem_ji10]);
		lctmp10 -= _OFFD_SP_clov[0][elem_ij10];
		lctmp10 -= adj(_OFFD_SP_clov[1][elem_ji10]);
		lctmp10 += _OFFD_SP_clov[1][elem_ij10];

		_SP_B[j10][i10] = lctmp10;

		elem_ij10 = elem_ij10 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;
    
  case 12:
    /*# gamma(  34)   i  0  0  0            # ( 1100 )  --> 12 */
    /*#               0 -i  0  0 */
    /*#               0  0 -i  0 */
    /*#               0  0  0  i */
    /*# From diagonal part */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp12;
	InnerReal lr_zero12;
	InnerReal lrtmp12;
	Complex sqm12;
	int i12;
	int j12;
	int elem_ij12;
	int elem_ijb12;
  
                          
	lr_zero12 = 0;
	FILL(sqm12, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i12 = 0; i12 < Nc; ++i12)
	    {
	      lrtmp12 = _DIAG_SP_clov[0][i12];
	      lrtmp12 -= _DIAG_SP_clov[0][i12+Nc];
	      lrtmp12 -= _DIAG_SP_clov[1][i12];
	      lrtmp12 += _DIAG_SP_clov[1][i12+Nc];
	      _SP_B[i12][i12] = cmplx(lr_zero12,lrtmp12);
	    }
    
	    /*# From lower triangular portion */
	    elem_ij12 = 0;
	    for(i12 = 1; i12 < Nc; ++i12)
	    {
	      elem_ijb12 = (i12+Nc)*(i12+Nc-1)/2 + Nc;
	      for(j12 = 0; j12 < i12; ++j12)
	      {
		lctmp12 = _OFFD_SP_clov[0][elem_ij12];
		lctmp12 -= _OFFD_SP_clov[0][elem_ijb12];
		lctmp12 -= _OFFD_SP_clov[1][elem_ij12];
		lctmp12 += _OFFD_SP_clov[1][elem_ijb12];
	
		_SP_B[j12][i12] = lctmp12 * sqm12;
		_SP_B[i12][j12] = adj(lctmp12) * sqm12;
	
		elem_ij12 = elem_ij12 + 1;
		elem_ijb12 = elem_ijb12 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;
    
  default:
    B = 0;
  }
  

  END_CODE("subroutine");;
}
#endif
