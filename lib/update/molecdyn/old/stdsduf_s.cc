/* $Id: stdsduf_s.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/* This routine is specific to staggered fermions! */

/* StDsDuf -- computes the derivative of the fermionic action respect the  */
/*            link field */
/* u -- gauge field ( Read ) */
/*         |  dS      dS_f */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU */
/* psi -- [1./(M_dag*M)]*chi_  ( read ) */

void StDsDuf(multi1d<LatticeColorMatrix>& ds_u,
	     const multi1d<LatticeColorMatrix>& u,
	     const LatticeFermion& psi)
{
  LatticeFermion u_psi;
  LatticeFermion rho;
  LatticeFermion u_rho;

  LatticeFermion tmp_1;
  Real dummy;
  int mu;
  int cb;

  START_CODE("subroutine");;

  /* rho = Dslash(1<-0) * psi */
  dslash (u, psi, rho, PLUS, 0);

  /* rho = (KappaMD^2)*rho = (KappaMC^2)*Dslash*psi */
  dummy = KappaMD*KappaMD;
  rho = rho * dummy;

  for(mu = 0;mu  < ( Nd); ++mu )
  {
    cb = 0;

    /* tmp_1(x) = rho(x+mu) */
    tmp_1(rb[cb]) = shift(rho, FORWARD, mu);

    /* u_rho(x) = u(x,mu)*rho(x+mu) = u(x,mu)*tmp_1(x) */
    /* Note the KS phase factors are already included in the U's! */
    u_rho = u[mu][cb] * tmp_1;
    
    /* ds_u(x,mu) = - u_rho(x) * psi_dag(x,mu) */
    ds_u[mu][cb] -= u_rho * adj(psi);
        
    cb = 1;
    
    /* tmp_1(x) = psi(x+mu) */
    tmp_1(rb[cb]) = shift(psi, FORWARD, mu);
    
    /* u_psi(x) = u(x,mu)*psi(x+mu) = u(x,mu)*tmp_1(x) */
    /* Note the KS phase factors are already included in the U's! */
    u_psi = u[mu][cb] * tmp_1;
    
    /* ds_u(x,mu) = u_psi(x) * rho_dag(x,mu) */
    ds_u[mu][cb] += u_psi * adj(rho);
  }
      
  END_CODE("subroutine");;
}
