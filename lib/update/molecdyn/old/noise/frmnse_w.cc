// $Id: frmnse_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO FIT TO USE VIRTUAL FUNCS"

#include "chromabase.h"


/* This routine is specific to Wilson fermions! */

//! Calculates the Pseudofermion fields and the pseudofermionic action
/*!
 * \ingroup molecdyn
 *
 */
/* u           -- gauge field ( Read ) */
/* chi         -- pseudofermion field. ( Write ) */
/* psi         -- pseudofermion field. ( Write ) */
/* chi_norm    -- | chi | ( Write ) */
/* n_zero      -- Number and type of zero modes (for Overlap) (Read) */
/* Ncb         -- Number of checkerboards (Read) */
/* Npf         -- Number of pseudofermion fields (Read) */
/* n_congrd    -- Number of CG iteration ( Write ) */

/* eta = Gaussian random noise */
/* chi = M_dag * eta */
/* psi = (M_dag*M)^(-1) * chi */

void FrmNse(multi1d<LatticeColorMatrix>& u,
	    multi1d<LatticeFermion>& chi,
	    multi1d<LatticeFermion>& psi,
	    int Npf,
	    int& n_congrd;
	    int n_zero)
{
  START_CODE();
  
  int n_count;
  n_congrd = 0;

  for(int i = 0; i < Npf; ++i)
  {
    switch (FermAct)
    {
    case WILSON:
    case PARITY_BREAKING_WILSON:
    case CLOVER:
      /* The noise is generated the same for these fermions */
      WlFrmNse(u, chi[i][0], psi[i][0], n_count);
      break;
    
    case OVERLAP_POLE:
      OvFrmNse(u, chi[i][0], psi[i][0], n_zero, n_count);
      break;

    default:
      QDP_error_exit("Unsupported fermion action", FermAct);
    }

    n_congrd += n_count;
  }

  END_CODE();
}
