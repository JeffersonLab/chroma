/* $Id: hywrmup_orig.cc,v 1.1 2003-12-31 23:48:22 edwards Exp $ ($Date: 2003-12-31 23:48:22 $) */

/* HYWRMUP */

/* Routine for warming up the gauge field configuration configuration  */
/* for the hybrid (monte carlo) algorithm. */

/* Modifications: */
/*   (ROSSI/ADK/RGE) Took the section out of szinwl to do this procedure. */
/* Sat May  5 20:25:55 EDT 1990 (ROSSI/ADK/RGE) */
/* Mon Jul 16 21:27:05 EDT 1990 (RGE) Changed implicit statement. */
/* Tue Oct 15 23:28:34 EDT 1991 (ADK/RGE) Changed traj. update to call to */
/*                              HybTrj. */
/* Wed Jan  8 17:42:38 EST 1992 (RGE) Made RsdCGMD an input parameter to invert  */
include(types.mh)

SUBROUTINE(HyWrmUp, u, WarmUp, Ncb, Npf)

/* Start declarations */
/* Argument declarations */

multi1d<LatticeColorMatrix> u(Nd);	/* New gauge field cfg. */
int WarmUp;                /* Warm-up (equilibration) trajectories */
int Ncb;                   /* Number of checkerboards */
int Npf;                   /* Number of Pseudo-Fermions */

{ /* Local variables */

  include(COMMON_DECLARATIONS)
  multi1d<LatticeColorMatrix> p_mom(Nd);  /* Fictitious momenta  */
  multi2d<LatticeFermion> chi(Ncb, Npf);			      /* Pseudofermion field */
  multi2d<LatticeFermion> psi(Ncb, Npf);			      /* Inv[Dag(M).M] Chi   */
  multi1d<LatticeFermion> phi(Ncb);

  int oldalg;		/* old AlgLPStp (temporary) */
  int n_loop;		/* Loop counter over trajectories */
  int cg_iter;		/* CG iterations for current trajectory */
  int n_count;		/* CG iterations performed by subroutine */
  int n_acc;		/* Number of MC acceptances */
  int i;                   /* Dummy - counter */       
  int ichiral;             /* Chirality of zero mode */
  int cb;

  Double ke;                     /* Kinetic energy */
  Double old_ke;                 /* ... of old CFG */
  Double pe;                     /* Bosonic action (potential energy) */
  Double old_pe;                 /* ... of old CFG */
  Double fe;                     /* Fermionic action (fermion energy) */
  Double old_fe;                 /* ... of old CFG */
  multi1d<Double> chi_norm(Npf);

  Double temp;
  Double temp_2;

  /* Include any hooks in the start macro */
  START_CODE("subroutine");;			/*# May include declarations and/or code. */

  /* Start executable code */
  /* Standard temporaries */
        psi = 0;

  /* For warmup we refresh momenta and fermion noise for every trajectory */
  RefNextTrj = YES;

  /* Perform warm-up trajectories */
  for(n_loop=1; n_loop <= WarmUp; ++n_loop)
  {
    
    if ( FermAct == OVERLAP_POLE && FermiP == YES)
      {
	if ( NOperEigDim > 0 )
	  CalcEigVec (u, NO);

	ZeroSector (u, 1, ichiral, Ncb);
	ichiral = - ichiral;
      }
    /* Generate momenta from Gaussian heatbath */
    Refrsh (p_mom);

    /* Generate pseudofermions from heatbath */
    if ( FermiP == YES)
    { 
      if ( PolyEvolP == YES )
      {
	for(i = 0; i < Npf; ++i)
	{
	  PolyFrmNse (u, chi[i][0], Ncb, n_count);
	  cg_iter = cg_iter + n_count;
	}
	chi_norm = 0;
      }
      else if ( RatEvolP == YES )
      {
		for(i = 0; i < Npf; ++i)
	  {
	    for (cb=0; cb<Ncb; ++cb) phi[cb] = chi[i][cb];
	    RatFrmNse (u, phi, Ncb, n_count);
	    for (cb=0; cb<Ncb; ++cb) chi[i][cb] = phi[cb];	    
	    cg_iter = cg_iter + n_count;
	    /*    temp = norm2(chi[0][0]);
		  PRINTF("%e",temp);
		  FLUSH_WRITE_NAMELIST(stdout);
	    	    eta_norm = 0;
		    for (cb=0; cb<Ncb; ++cb) {
		    eta_norm += norm2(psi[cb]);
		    }
		    eta_norm = sqrt(eta_norm);
		    printf("%e",eta_norm);*/
	  }
	chi_norm = 0;
	      }
      else
      {
	FrmNse (u, chi, psi, chi_norm, ichiral, Ncb, Npf, n_count);
	cg_iter = cg_iter + n_count;
      }
    }
    else
    {
      if ( FermiP == NO ) chi_norm = 0;
    }

    /* Measure the contributions to the MD Hamiltonian */
    MesE (u, p_mom, chi, psi, BetaMD, KappaMD, old_ke, old_pe, old_fe, Ncb, Npf, n_count);
    cg_iter = cg_iter + n_count;
        

    /* One complete HMD trajectory */
    ke = old_ke;
    pe = old_pe;
    fe = old_fe;
    
    oldalg   = AlgLPStp;
    AlgLPStp = NO;
    HybTrj (u, p_mom, chi, chi_norm, psi, n_count, ke, pe, fe, Ncb, Npf);
    cg_iter = cg_iter + n_count;
    AlgLPStp = oldalg;
    
    old_ke = ke;
    old_pe = pe;
    old_fe = fe;

    /* Do not perform last half-step, because we do not need final momenta during */
    /* warm-up as no Metropolis step is used */
      }				/* end loop for warm-up trajectories */

        /* Close out any other code */
  END_CODE("subroutine");;
}

