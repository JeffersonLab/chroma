/* $Id: dsduf_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/* This routine is specific to Wilson fermions! */

/* dsduf -- computes the derivative of the fermionic action respect the  */
/*          link field */
/* u -- gauge field ( Read ) */

/*         |  dS      dS_f */
/* ds_u -- | ----   + -----   ( Modify ) */
/*         |  dU       dU */

/* chi      -- pseudofermion field. ( Read ) */
/* chi_norm -- | chi | ( Read ) */
/* psi -- [1./(M_dag*M)]*chi_  ( read ) */
/* Ncb   --  Number of checkerboards. Ignored for now     (Read) */
/* Npf   -- number of pseudofermions  ( Read ) */
/* cg_count -- number of CG iterations in force routine ( Modify ) */

include(types.mh)

SUBROUTINE(dsduf, u, ds_u, chi, chi_norm, psi, Ncb, Npf, cg_count)

multi1d<LatticeColorMatrix> u(Nd);
multi1d<LatticeColorMatrix> ds_u(Nd);
multi2d<LatticeFermion> psi(Ncb, Npf);
multi2d<LatticeFermion> chi(Ncb, Npf);
multi1d<Double> chi_norm(Npf);
int Ncb;
int Npf;
int cg_count;
{ /* Local variables */
  include(COMMON_DECLARATIONS)

  int cb;
  int i;

  START_CODE("subroutine");;
  
  cg_count = 0;

  for(i = 0; i < Npf; ++i)
  {
    switch (FermAct)
    {
    case WILSON:
      /* Do the usual Wilson fermion dS_f/dU */
      WlDsDuf (u, ds_u, psi[i][0]);
      break;
    
    case PARITY_BREAKING_WILSON:
      /* Do the Wilson fermion dS_f/dU with parity breaking term */
      WlhDsDuf (u, ds_u, psi[i][0]);
      break;
    
    case CLOVER:
      /* Do the clover fermion dS_f/dU */
      ClovDsDuf (u, ds_u, psi[i][0]);
      break;
      
    case OVERLAP_POLE:
      /* Do the overlap fermion dS_f/dU */
      OvDsDuf (u, ds_u, psi[i][0], Ncb);
      break;
    
    default:
      QDP_error_exit("Unknown fermion action", FermAct);
    }
  }

  END_CODE("subroutine");;
}
