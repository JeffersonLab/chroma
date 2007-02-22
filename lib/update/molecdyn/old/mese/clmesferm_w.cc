// $Id: clmesferm_w.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $

#error "NOT FULLY CONVERTED - MAJOR CHANGES OR SUPPORT FOR FIELD-STRENGTH F NEEDED. HERE FOR IDEAS TO HELP WITH OVERALL CLASS STRUCTURE"

#include "chromabase.h"


/* This routine is specific to Wilson fermions! */

/* MESFERM -- Measures the pseudofermion contribution to the Hamiltonian. */
/*            THIS CODE USES THE SOLUTION PSI PASSED TO IT. IT DOES NOT */
/*            RECOMPUTE THE CG SOLUTION.  */

/* u      -- gauge field (Read) */
/* chi    -- Gaussian random noise for pseudo-fermion source (Read) */
/* psi    -- (M_dagM)^(-1) * chi (Read) */
/* Kappa  -- Kappa (Read) */
/* w_ferm -- Pseudo fermion energy (Write) */
/* n_congrd  --  number of CG iterations (Write)                          */
/* NOTES: */
/* w_ferm = chi_dag * psi = chi_*(M_dagM)^(-1)*chi = |(1/M_dag)*chi|^2 */

void ClMesFerm(const multi1d<LatticeColorMatrix>& u,
	       const LatticeFermion& chi,
	       const LatticeFermion& psi,
	       const Real& Kappa,
	       Double& w_ferm,
	       const Subset& sub)
{   
  LATTICE_TRIANG(clov);
  LATTICE_TRIANG(invclov);
  LATTICE_FIELD_STRENGTH(f);
  LatticeReal trace_aux;
  Double logdet;
  Real Kappa_cl;
  
  START_CODE();
  
  if ( FermiP == YES )
  {
    /* chi_bar*(1/(M_dag*M))*chi = chi_bar*psi */
    w_ferm = innerProduct(chi,psi,sub) / Double(Layout::vol()*Nc*Ns);

    /* compute also contribution from the "preconditioning matrix" */
            
    /*+ */
    /* The Kappas used in dslash and clover mass term */
    /* Kappa_ds = Kappa_dslash = Kappa */
    /* Kappa_cl = Kappa_clovms = ClovCoeff * Kappa / u0^3 */
    /*- */
    Kappa_cl = ClovCoeff * Kappa / (u0*u0*u0);
    if (AnisoP == YES)
      QDP_error_exit("anisotropy not supported");

    /* Calculate F(mu,nu) */
    MesField(u, f);

    /* Make the 'clover term', and compute the log(det) */
    makclov(f, clov, Kappa_cl, Kappa_cl, 0);
    chlclovms(clov, invclov, YES, logdet);

    /* Normalize by number of pseudo-fermionic degree of freedom */
    logdet /= Double(Layout::vol()*Nc*Ns);
    w_ferm -= 2*logdet;   // HMM, THE 2 IS THE NUMBER OF FERMIONS??? DON'T RECALL...
  }
  else
  {
    w_ferm = 0;
  }
  
  END_CODE();
}
