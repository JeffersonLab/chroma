/* $Id: dsduf_s.cc,v 2.0 2005-09-25 21:04:42 edwards Exp $ ($Date: 2005-09-25 21:04:42 $) */

/* This routine is specific to staggered fermions! */

/* dsduf -- computes the derivative of the fermionic action respect the  */
/*          link field */
/* u -- gauge field ( Read ) */
/*         |  dS      dS_f */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU */
/* psi -- [1./(M_dag*M)]*chi_  ( read ) */
/* Ncb   --  Number of checkerboards. Ignored for now     (Read) */
/* Npf   -- number of pseudofermions  ( Read ) */
/* cg_count -- number of CG iterations in force routine ( Modify ) */

void dsduf(multi1d<LatticeColorMatrix>& ds_u,
	   const multi1d<LatticeColorMatrix>& u,
	   const multi1d<LatticeFermion>& psi,
	   const multi1d<LatticeFermion>& chi,
	   int Npf,
	   int cg_count)
{
  LatticeFermion tmp;
  int i;

  START_CODE("subroutine");;

  if (FermAct != STAGGERED)
    QDP_error_exit("unsupported fermion action", FermAct);

  if (Ncb != 1)
    QDP_error_exit("expect Ncb = 1", Ncb);

  cg_count = 0;

  for(i = 0; i < Npf; ++i)
  {
    StDsDuf (u, ds_u, psi[i][0]);
  }
  
  END_CODE("subroutine");;
}
