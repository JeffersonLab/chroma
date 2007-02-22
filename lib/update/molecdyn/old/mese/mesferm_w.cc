// $Id: mesferm_w.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $

#error "NOT FULLY CONVERTED"

#include "chromabase.h"


/* This routine is specific to Wilson fermions! */

//! MESFERM -- Measures the pseudofermion contribution to the Hamiltonian.
/*!
 * \ingroup molecdyn
 *
 */

/*            THIS CODE USES THE SOLUTION PSI PASSED TO IT. IT DOES NOT */
/*            RECOMPUTE THE CG SOLUTION.  */

/* u      -- gauge field (Read) */
/* chi    -- Gaussian random noise for pseudo-fermion source (Read) */
/* psi    -- (M_dagM)^(-1) * chi (Read) */
/* Kappa  -- Kappa (Read) */
/* w_ferm -- Pseudo fermion energy (Write) */
/* n_congrd  --  number of CG iterations (Write)                          */
/* Npf    -- number of pseudofermion fields  ( Read ) */
/* NOTES: */
/* w_ferm = chi_dag * psi = chi_*(M_dagM)^(-1)*chi = |(1/M_dag)*chi|^2 */

void MesFerm(const multi1d<LatticeColorMatrix>& u,
	     const multi1d<LatticeFermion>& chi,
	     const multi1d<LatticeFermion>& psi,
	     const Real& Kappa,
	     Double& fe,
	     int& n_congrd,
	     int Npf,
	     const Subset& sub)
{
  START_CODE();
  
  n_congrd = 0;
  fe = 0;

  if ( FermiP == YES )
  {
    Double w_ferm;

    for(int i = 0; i < Npf; ++i)
    {
      switch (FermAct)
      {
      case WILSON:
      case PARITY_BREAKING_WILSON:
      case OVERLAP_POLE:
	/* The noise is generated the same for these fermions */
	WlMesFerm(u, chi[i], psi[i], w_ferm, sub);
	break;
    
      case CLOVER:
	ClMesFerm(u, chi[i], psi[i], Kappa, w_ferm, sub);
	break;

      default:
	QDP_error_exit("Unsupported fermion action", FermAct);
      }

      fe += w_ferm;
    }
  }
  
  END_CODE();
}
