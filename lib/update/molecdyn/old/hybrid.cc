// $Id: hybrid.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $
// HYBRID

#error "NOT FULLY CONVERTED"

#include "chromabase.h"


//! Routine for doing the hybrid (monte carlo) algorithm.
/*!
 * \ingroup molecdyn
 *
 */
/* u        -- gauge field ( Modify ) */
/* Algorithm -- flag for Monte Carlo step ( Read ) */
/* MonitordH -- Monitor dH along a trajectory ( Read ) */
/* t_srce   -- Cartesian source coordinates for spectroscopy ( Read ) */
/* j_decay  -- Cartesian source coordinates for spectroscopy ( Read ) */
/* MesonP   -- flag for spectroscopy ( Modify ) */
/*             (0) no measurement. */
/*             (1) for measurement iff configuration is *really* accepted */
/*                 (the first configuration of a batch is always accepted). */
/*             (-1) will force measurement. Is reset to "1" at end of  */
/*                 measurement. */
/* GlueP  -- flag for glueball measurements ( Read ) */
/* WilsonP  -- flag for Wilson loop measurements ( Read ) */
/* SchrPcacP -- flag for Schroedinger PCAC measurements ( Read ) */
/* MesTrj   -- number of trajectories in mod for measurements ( Read ) */
/* SpecTrj  -- number of trajectories in mod for spectroscopy ( Read ) */
/* NumCG    -- number of CG iterations performed ( Write ) */
/* TotalCG  -- total number of CG iterations performed ( Modify ) */
/* NumTrj   -- number of trajectories to perform ( Read )  */
/* TotalTrj -- total number of trajectories to perform ( Read )  */
/* RefMomTrj -- trajectories per momentum refreshment ( Read )  */
/* RefFnoiseTrj -- trajectories per fermion noise refreshment ( Read )  */
/* spec_acc -- True if new cfg accepted since spectroscopy meas. ( Modify ) */
/* BlkAccu  -- blocking accuracy in glueball code ( Read ) */
/* BlkMax   -- maximum number of SU(3) trace maximimations allowed ( Read ) */
/* GFAccu   -- gauge fixing accuracy ( Read ) */
/* GFMax    -- maximum number of gauge fixing relaxation steps ( Read ) */
/* OrDo     -- flag for over-relaxation ( Read ) */
/* OrPara   -- over-relaxation parameter ( Read ) */
/* TopolP   -- flag for topology measurements ( Read )# */
/* TopAccu  -- accuracy for convergence of topological charge ( Read ) */
/* ActAccu  -- accuracy for convergence of action ratio ( Read ) */
/* NumTop   -- number of topological charge measurements ( Read ) */
/* NumCool  -- number of cooling sweeps per topological charge measurement ( Read ) */
/* SmwlsP   -- Do "smeared" Wilson loop measurements? */
/* numb_sm  -- Number of smearing levels for "smeared" loops */
/* fact_sm  -- 'Smearing' factor for "smeared" loops */
/* Npf      -- number of pseudofermions  ( Read ) */

void hybrid(u, 
	     t_srce, j_decay, MesonP, GlueP, WilsonP, TopolP, SchrPcacP, spec_acc,
	     SmwlsP, numb_sm, fact_sm,
	     Hybrid_params, NumTrj, TotalTrj,
	     NumCG, TotalCG, 
	     BlkMax, BlkAccu, GFAccu, GFMax,
	     OrDo, OrPara,
	     TopAccu, ActAccu, NumTop, NumCool, Npf)

/*  PrtMes, MesTrj, SpecTrj, RefMomTrj, RefFnoiseTrj  */

/* Start declarations */
/* Argument declarations */

  multi1d<LatticeColorMatrix> u(Nd);	/* New gauge field cfg. */

multi1d<int> t_srce(Nd);		/* Coordinates of the quark propagator source. */
int j_decay;		/* Direction to measure propagators. */
int MesonP;		/* Measure meson propagators on-line? */
int GlueP;			/* Measure glueballs on-line? */
int WilsonP;		/* Measure Wilson loops on-line? */
int NumCG;			/* Number of CG iterations for this run */
int TotalCG;		/* Total number of CG iterations */
int NumTrj;		/* Number of trajectories for this run */
int TotalTrj;		/* Total number of trajectories */
int spec_acc;		/* Acceptance flag for spectroscopy */

multi1d<int> Hybrid_params(5);       /* Really stupid way to pack some variables to avoid 
				 * too many arguments for hybrid */ 

int SmwlsP;		/* Do "smeared" Wilson loop measurements? */
int numb_sm;		/* Number of smearing levels for "smeared" loops */
Real fact_sm;			/* 'Smearing' factor for "smeared" loops */

int TopolP;		/* Measure topology on-line? */
int SchrPcacP;
Real TopAccu;                  /* Topological charge measurement accuracy */
Real ActAccu;                  /* Action measurement accuracy */
int NumTop;		/* Frequency of charge measurement */
int NumCool;		/* Total number of cooling sweeps */

Real BlkAccu;		        /* Accuracy of blocking in fuzglue. */
int BlkMax;		/* Maximum number of proj. iterations in fuzglue. */
Real GFAccu;			/* Accuracy of gauge fixing */
int GFMax;			/* Maximum number of gauge-fixing relaxations */
int OrDo;			/* Do over-relaxation? */
Real OrPara;			/* Value of over-relaxation parameter (OrDo=YES) */
int Ncb;                   /* Number of checkerboards */
int Npf;                   /* Number of pseudofermion fields */

{
  /* Variables on input but packed into Hybrid_params */
  int PrtMes;		/* Measurements per printout to RESULT file */
  int MesTrj;		/* Trajectories per measurement */
  int SpecTrj;		/* Trajectories per spectroscopy */
  int RefMomTrj;		/* Trajectories per momentum refreshment */
  int RefFnoiseTrj;	/* Trajectories per fermion noise refreshment */

  /* Local variables */ 
  multi1d<LatticeColorMatrix> u_tmp(Nd);	/* Old gauge field cfg.  */
  multi1d<LatticeColorMatrix> u_smr(Nd);	/* Smeared cfg. */
  multi1d<LatticeColorMatrix> p_mom(Nd);	/* Fictitious momenta */
  multi2d<LatticeFermion> chi(Ncb, Npf);				/* Pseudofermion field */
  multi2d<LatticeFermion> psi(Ncb, Npf);				/* Inv[Dag(M).M] Chi */
  
  int mu;
  
  multi1d<DComplex> pollp(Nd);		/* Polyakov loop */
  Double w_plaq;			/* Whole plaquette */
  Double s_plaq;			/* Spatial plaquette */
  Double t_plaq;			/* Temporal (thermal) plaquette */
  Double link;			/* Trace of link (not gauge invariant) */
  
  int ReverseChk;	        /* Flag indicating reversibility check */
  int Reverse_loop;	        /* Do loop counter for reversibility check */
  int Reverse_num;	        /* Number of loops for reversibility check */
  int TrjNumber;	        /* Junk variable indicating the trajectory */
  int n_loop;		/* Loop counter over trajectories */
  int cg_iter;		/* CG iterations for current trajectory */
  int n_count;		/* CG iterations performed by subroutine */
  int n_acc;			/* Number of MC acceptances */
  int conf_acc;		/* Was new CFG accepted? */
  int NProbViLa;           /* # of probability estimate violations > 1 */
  int NProbViSm;           /* # of probability estimate violations < 0 */
  int NTermsExp;           /* # of nontrivial terms to estimate Exp    */
  int NTermsLog;           /* # of nontrivial terms to estimate Log    */
  int ncg_had;		/* Number of CG iterations in hadron */
  int nrl_gf;		/* Number of relaxations in gfix */
  int length;		/* Length along direction of decay */
  int half_length;		/* Half-length along direction of decay */
  int t0;			/* Coordinate of slice in direction of decay */
  multi1d<int> t_srce_tmp(Nd);		/* Temp for hadron source coordinates */
  int numKappa;		/* Number of Kappa's passed to hadron */
  int bc_spec;		/* Boundary condition for spectroscopy */
  multi1d<Real> Kappa(numKappa);		/* Array of Kappa's passed to hadron */
  multi1d<Real> ClCoeff(numKappa);	/* Array of clover coeffs */
  
  Double ke;			/* Kinetic energy */
  Double old_ke;			/* ... of old CFG */
  Double rvse_ke;		/* ... for reversibility checks */
  Double guide_ke;		/* ... for guidance checks */
  Double pe;			/* Bosonic action (potential energy) */
  Double old_pe;			/* ... of old CFG */
  Double rvse_pe;		/* ... for reversibility checks */
  Double guide_pe;		/* ... for guidance checks */
  Double fe;			/* Fermionic action (fermion energy) */
  Double old_fe;			/* ... of old CFG */
  Double rvse_fe;		/* ... for reversibility checks */
  Double guide_fe;		/* ... for guidance checks */
  multi1d<Double> chi_norm(Npf);
  Double pbp_st;		/* Stochastic estimator of chiral condensate */
  Real DelH;			/* Change in total energy */
  Real DelHB;			/* Change in bosonic energy */
  Real DelKe;			/* Change in kinetic energy */
  Real DelPe;			/* Change in potential energy */
  Real DelFe;			/* Change in fermionic energy */
  Real rand;			/* Uniform pseudo-random number in [0,1] */
  Real SumAccProb;		/* Accumulator for min(1,exp(-dH)) */
  Real AccProb;   		/* Acceptance probability = min(1,exp(-dH)) */
  Real DetRat;                 /* Estimate of the determinant ratio   */
  Real lam_lo;                 /* Lowest eigenvalue (for MesEV call)  */
  Real lam_hi;                 /* Highest eigenvalue (for MesEV call) */
  
  int ichiral;             /* Chirality of zero mode */

  multi1d<LatticeFermion> phi(Ncb);
  Double rdummy;                 /* Scratch variable */
  Double idummy;                 /* Scratch variable */
  int cb;                    /* Scratch variable */
  int i;                     /* Scratch variable */
  Real ftmp1;                    /* Scratch variable */
  Real ftmp2;                    /* Scratch variable */
  int itmp1;                 /* Scratch variable */
  
  /* Include any hooks in the start macro */
  START_CODE("subroutine");;			/*# May include declarations and/or code. */
  
  /* Unpack variables */
  PrtMes = Hybrid_params[0];
  MesTrj = Hybrid_params[1];
  SpecTrj = Hybrid_params[2];
  RefMomTrj = Hybrid_params[3];
  RefFnoiseTrj = Hybrid_params[4];

  /* Start executable code */
  /*# The starting initial cfg acceptance is *NOT* used as a flag here */
  /* for measurements. Measurements are always done depending on MesTrj. */
  /* Spectroscopy is done based on spec_acc, MesonP and SpecTrj. When MesonP is */
  /* first turned on (1), the spectroscopy is forced to be measured when */
  /* mod(TotalTrj,SpecTrj) = 0. Otherwise, it is only done on acceptance,  */
  /* and spec_acc == YES, on MesonP == 1, and when mod(TotrlTrj,SpecTrj) = 0. */
  /* NOTE: spec_acc is a latched on when any configuration was accepted over */
  /* SpecTrj trajectories. Spec_acc is turned off after the measurements. */
  SumAccProb = 0;
  conf_acc = NO;
  
  /* Compile time option to handle reversibility check */
  ifdef(`HYBRID_REVERSIBILITY_CHECK',`ReverseChk = YES; MonitordH = YES',dnl
  `ReverseChk = NO')`';

  /* On-line statistics */
  pollp = 0;

  n_acc = 0;
  w_plaq = 0;
  s_plaq = 0;
  t_plaq = 0;
  link = 0;
  DelH = 0;

  /* Initialize the counters for stochastic probability estimates in R-mode */
  NProbViLa = 0;
  NProbViSm = 0;
  NTermsExp = 0;
  NTermsLog = 0;
  
  /* Perform Hybrid Monte Carlo algorithm */

        
  chi = 0;
  psi = 0;

  /* Do NumTrj HMC trajectories */
  for(n_loop = 1; n_loop <= NumTrj; ++n_loop )
  {
    RoundP = NO;
    cg_iter = 0;		/* CG iterations for current trajectory */

    u_tmp = u;	/* keep old CFG for Metropolis rejection */

    TotalTrj = TotalTrj + 1;    /* Trajectory bookkeeping.     */
    TrjNumber = TotalTrj;	/* Want a different name here. */

    push(nml_out,"NewTrajectory");
    Write(nml_out, TrjNumber);
    pop(nml_out);

    if ( FermAct == OVERLAP_POLE && FermiP == YES )
    {
      if ( ( n_loop == 1 ) || 
	   (( ((TotalTrj-1) % MesTrj)  != 0 ) &&
	    ( ((TotalTrj-1) % SpecTrj) != 0 ))) 
      {
	if ( NOperEigDim > 0 )
	  CalcEigVec (u_tmp, NO);

	ZeroSector (u_tmp, 5, ichiral, Ncb);
	ichiral = - ichiral;
      }
    }
    else if ( FermAct == TRUNC_OVERLAP && FermiP == YES )
    {
      if ( ( n_loop == 1 ) || 
	   (( ((TotalTrj-1) % MesTrj)  != 0 ) &&
	    ( ((TotalTrj-1) % SpecTrj) != 0 ))) 
      {
	if ( NOperEigDim > 0 )
	  CalcEigVec (u_tmp, NO);
      }
    }

    /* Generate momenta from Gaussian heatbath */
    if ( (TotalTrj % RefMomTrj) == 0 )
      Refrsh (p_mom);

    /* Generate pseudofermions from heatbath */
    if ( FermiP == YES && (TotalTrj % RefFnoiseTrj) == 0 &&
	 Algorithm != R)
    { 
      if ( PolyEvolP == YES )
      {
	for(i = 0; i < Npf; ++i)
	{
	  PolyFrmNse (u_tmp, chi[i][0], Ncb, n_count);
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
	}

	chi_norm = 0;
      }
      else
      {
	FrmNse (u_tmp, chi, psi, chi_norm, ichiral, Ncb, Npf, n_count);
	cg_iter = cg_iter + n_count;
      }
    }
    else
    {
      if ( FermiP == NO ) chi_norm = 0;
    }

    /* Will P or chi be refreshed for next trajectory?         */
    if ( (TotalTrj % RefFnoiseTrj) == 0 ||                  \
         (TotalTrj % RefMomTrj) == 0 )
      RefNextTrj = YES;
    else
      RefNextTrj = NO;      	               
       
    /* Start to calculate initial Hamiltonian (we could be more economical  */
    /*# and not recompute the plaquettes, but this is a negligible saving)  */
    /* Measure the contributions to the MC Hamiltonian. */

    MesE (u_tmp, p_mom, chi, psi, BetaMC, KappaMC, old_ke, old_pe, old_fe, Ncb, Npf, n_count);
    cg_iter = cg_iter + n_count;    

    if ( KappaMC != KappaMD || BetaMC != BetaMD || RsdCGMC != RsdCGMD )
    {
      if ( KappaMC != KappaMD || RsdCGMC != RsdCGMD )
      {
	if ( ( FermiP == YES ) && ( PolyEvolP == NO ) && ( RatEvolP == NO ) )
	{
	  /* psi = (1/(M_dag*M))*chi */
	  phfctr (u_tmp, FORWARD);	   /* ON */
	  psi = 0;
          for( i = 0; i < Npf; ++i )
	  {
            invert (u_tmp, chi[i][0], psi[i][0], chi_norm[i], KappaMD, RsdCGMD, Ncb, n_count);
	    cg_iter = cg_iter + n_count;
          }
	  phfctr (u_tmp, BACKWARD);   /* OFF */
	}
      }

      MesE (u_tmp, p_mom, chi, psi, BetaMD, KappaMD, ke, pe, fe, Ncb, Npf, n_count);
      cg_iter = cg_iter + n_count;    

      guide_ke = ke;
      guide_pe = pe;
      guide_fe = fe;
    }
    else
    {
      ke = old_ke;
      pe = old_pe;
      fe = old_fe;
    }

    /* Perform a trajectory specified by the flag Algorithm */
    if ( ReverseChk == YES )
    {
      rvse_ke = ke;
      rvse_pe = pe;
      rvse_fe = fe;
      push(nml_out,"Initial_Reversibilty_Check");
      Write(nml_out, rvse_ke);
      Write(nml_out, rvse_pe);
      Write(nml_out, rvse_fe);
      pop(nml_out);
      Reverse_num = 2;
    }
    else
      Reverse_num = 1;


    for(Reverse_loop = 1;Reverse_loop  <= ( Reverse_num); ++Reverse_loop )
    {
      RoundP = YES;

      if ( Algorithm == HMCC || Algorithm == HMDC )
      {
	/*     Perform a Campostrini style trajectory */
	CampTrj (u_tmp, p_mom, chi, chi_norm, psi, n_count, ke, pe, fe, Ncb, Npf);
	cg_iter = cg_iter + n_count;
      }
      else if ( Algorithm == R) {
	/*     Perform an R algorithm trajectory.*/
	RHybTrj (u_tmp, p_mom, chi, chi_norm, psi, n_count, ke, pe, fe, Ncb, Npf);
	cg_iter = cg_iter + n_count;
      } 
      else 
      {
	/*     This is a standard single forward step trajectory */
	HybTrj (u_tmp, p_mom, chi, chi_norm, psi, n_count, ke, pe, fe, Ncb, Npf);
	cg_iter = cg_iter + n_count;
      }

      /*+ */
      /* For reversibility check, flip the momenta and continue. This will */
      /* integrate in the reverse direction the equations of motion of the */
      /* current configuration *on the *next* time through the loop*.      */
      /*- */
      if ( ReverseChk == YES )
      {
	DelKe = ke - rvse_ke;
	DelPe = pe - rvse_pe;
	DelFe = fe - rvse_fe;
	DelH  = DelKe - DelPe + DelFe;
	push(nml_out,"Reversibilty_Check");
	Write(nml_out, Reverse_loop);
	Write(nml_out, DelH);
	Write(nml_out, DelKe);
	Write(nml_out, DelPe);
	Write(nml_out, DelFe);
	pop(nml_out);
	FPRINTF(trm_out,"Reversibility check at TrjNumber = %d\n", TrjNumber);

	FlipMom (p_mom);
      }

      RoundP = NO;
    }                     /* End of Reverse_loop */


    /* Calculate the final Hamiltonian if we are using Hybrid Monte Carlo */
    /* Note that we use the AlgLPStp flag here, distinguishing the cases  */
    /* with/without Accept/Reject step                                    */
    if ( AlgLPStp == YES || MonitordH == YES )
    {
      if ( MonitordH == NO || BetaMC != BetaMD || KappaMC != KappaMD ||
	   RsdCGMD != RsdCGMC )
      {
	if ( BetaMC != BetaMD || KappaMC != KappaMD || RsdCGMD != RsdCGMC )
	{
	  MesE (u_tmp, p_mom, chi, psi, BetaMD, KappaMD, ke, pe, fe, Ncb, Npf, n_count);
          cg_iter = cg_iter + n_count;

	  /* Change in guidance energies */
	  DelKe = ke - guide_ke;
	  DelPe = pe - guide_pe;
	  DelFe = fe - guide_fe;
	  DelH  = DelKe - DelPe + DelFe;
	  push(nml_out,"Guidance_energy_diffs");
	  Write(nml_out, DelH);
	  Write(nml_out, DelKe);
	  Write(nml_out, DelPe);
	  Write(nml_out, DelFe);
	  pop(nml_out);
	}

	if ( FermAct == OVERLAP_POLE && FermiP == YES )
	{
	  ZeroSector (u_tmp, 5, ichiral, Ncb);
	}

	if ( KappaMC != KappaMD || RsdCGMD != RsdCGMC )
	{
	  /* New fermion solution required: More accurate or with KappaMC */
	  if ( ( FermiP == YES ) && ( PolyEvolP == NO ) && ( RatEvolP == NO ) )
	  {
	    phfctr (u_tmp, FORWARD);   /* ON */
	    psi = 0;
            for( i = 0; i < Npf; ++i )
            {
	      invert (u_tmp, chi[i][0], psi[i][0], chi_norm[i], KappaMC, RsdCGMC, Ncb, n_count);
	      cg_iter = cg_iter + n_count;
            }
	    phfctr (u_tmp, BACKWARD);  /* OFF */
	  }
	}

	/* Measure the contributions to the MC Hamiltonian */

	MesE (u_tmp, p_mom, chi, psi, BetaMC, KappaMC, ke, pe, fe, Ncb, Npf, n_count);
        cg_iter = cg_iter + n_count;
      }

      /* Find the change in energy per d.o.f. during the current trajectory */
      /* N.B.!! Potential energy appears with a minus sign! */

      DelKe = ke - old_ke;
      DelPe = pe - old_pe;
      DelFe = fe - old_fe;
      DelH  = DelKe - DelPe + DelFe;


      FPRINTF(trm_out,"     Old PE =%15.8g  Old KE =%15.8g  Old FE =%15.8g\n",
	      old_pe, old_ke, old_fe);
      FPRINTF(trm_out,"         PE =%15.8g      KE =%15.8g      FE =%15.8g\n     dH =%15.8g\n",
	      pe, ke, fe, DelH);

    } /* end of calculation of dH for HMC */

    /* Apply Metropolis accept/reject step. Even if we are using HMC,    */
    /* treat the case where DelH < 0 specially to avoid overflow in      */
    /* computing the exponential of the total energy change.             */
    /* Note that we use stochastic estimate of accept/reject probability */
    /* in the R-algorithm mode and in PHMCN algorithm.                   */

    /* Note that we use the AlgLPStp flag here, distinguishing the cases  */
    /* with/without Accept/Reject step                                    */
    if ( AlgLPStp == NO )
    {
      /* This is the branch without Accept/Reject test */
      /* Accept the new CFG */
      ftmp1 = 1;
      SumAccProb += ftmp1;
      u = u_tmp;
      n_acc = n_acc + 1;
      conf_acc = YES;
      spec_acc = YES;
    }
    else
    {
      /* This is the branch with Accept/Reject test */  
      if ( R0algP == NO && Algorithm != PHMCN && Algorithm != RHMCN)
      { 
	/* This is the branch with standard exact accept/reject test */
	if ( DelH < WORD_VALUE(WORD_DelH,ZERO) )
	{
	  AccProb = 1;
	  rand = 0;
	}
	else
	{
	  if ( Nc == 1)
	  {
	    AccProb = - TO_REAL(vol*Nd)*DelH;
	  }
	  else
	  { 
	    AccProb = - TO_REAL(vol*Nd*(Nc*Nc-1))*DelH;
          }
	  AccProb = exp(AccProb);
	  random(rand);
	}
	SumAccProb += AccProb;
	if ( rand <= AccProb )
	{
	  /* Accept the new CFG */
	  u = u_tmp;
	  n_acc = n_acc + 1;
	  conf_acc = YES;
	  spec_acc = YES;
	}
	else
	{
	  /* Reject the new CFG and keep the old one */
	  pe = old_pe; ke = old_ke; fe = old_fe;
	  /* Flip the momenta if they won't be refreshed */
	  if ( RefNextTrj == NO )  FlipMom (p_mom);
	}
      }
      else
      {  
	/* This is the branch with noisy accept/reject test */
	/* The stochastic estimate of the determinant ratio */
	EstDetRat (u, u_tmp, DetRat, n_count, NTermsExp, NTermsLog, Ncb);
	cg_iter = cg_iter + n_count;

	/* Nonstochastic part of acceptance probability  */
	DelHB = DelKe - DelPe;
	if (R0algP == YES)
	{
	  if ( Nc == 1)
	  {
	    AccProb = - TO_REAL(vol*Nd)*DelHB;
	  }
	  else
	  {
	    AccProb = - TO_REAL(vol*Nd*(Nc*Nc-1))*DelHB;
          }
	}
	else
	{
	  if ( Nc == 1)
	  {
	    AccProb = - TO_REAL(vol*Nd)*DelH;
	  }
	  else
	  {
	    AccProb = - TO_REAL(vol*Nd*(Nc*Nc-1))*DelH;
          }
	}
	AccProb = exp(AccProb);

	/* The estimate of total probability */       
	AccProb = AccProb*DetRat;

	/* The probability for Kennedy-Kuti acceptance step  */
	if ( DelH > WORD_VALUE(WORD_DelHB,ZERO) )
	{
	  AccProb = LamMi + LamPl*AccProb;
	}
	else
	{
	  AccProb = LamPl + LamMi*AccProb;
	}

	push(nml_out,"Stochastic_Stuff");
	Write(nml_out, DetRat);
	Write(nml_out, AccProb);
	pop(nml_out);

	/* Handle the violations of acceptance probability estimates */ 
	if (AccProb < 0.0)
	{ AccProb = 0;
	NProbViSm = NProbViSm + 1;
        push(nml_out,"Stoch_Exception");
	Write(nml_out, NProbViSm);
	pop(nml_out); 
	}
	if (AccProb > 1.0)
	{ AccProb = 1;
	NProbViLa = NProbViLa + 1;
        push(nml_out,"Stoch_Exception");
	Write(nml_out, NProbViLa);
	pop(nml_out);
	}

	SumAccProb += AccProb;
	random(rand);
	if ( rand <= AccProb ) 
	{                                   /* Accept the new CFG */
	  u = u_tmp;
	  n_acc = n_acc + 1;
	  conf_acc = YES;
	  spec_acc = YES;
	}
	else
	{                                   /* Reject the new CFG */
	  pe = old_pe; ke = old_ke; fe = old_fe;
	  /* Flip the momenta if they won't be refreshed */
	  if ( RefNextTrj == NO )  FlipMom (p_mom);
	}
      }
    }
    /* end of Metropolis accept/reject step */
    push(nml_out,"Cfg_acceptance");
    Write(nml_out, conf_acc);
    pop(nml_out);
  
      

    if ( (FermAct == OVERLAP_POLE || FermAct == TRUNC_OVERLAP) && 
         FermiP == YES &&
	 (( (TotalTrj % MesTrj) == 0 ) ||
          ( (TotalTrj % SpecTrj) == 0 ))) 
    {
      if ( (MonitordH == NO && AlgLPStp == NO && RefNextTrj == YES)
	   || (AlgLPStp == YES && conf_acc == NO) )
      {
	if ( NOperEigDim > 0 )
	  CalcEigVec (u, NO);
      }
      if ( FermAct == OVERLAP_POLE )
      {
	ZeroSector (u, 5, ichiral, Ncb);
	ichiral = - ichiral;
      }
    }
    /* Carry out measurements every MesTrj trajectories (typically MesTrj=1) */
    if ( (TotalTrj % MesTrj) == 0 )
    {
      MesPlq (u, w_plaq, s_plaq, t_plaq, link);

      if ( FermiP == YES )
      {

	if ( FermAct == PARITY_BREAKING_WILSON )
        {
          MesPbg5p(u, 1, pbp_st, n_count);
          cg_iter += n_count;
        }
        else
        {
	  MesPbp(u, pbp_st, n_count, ichiral);
	  cg_iter += n_count;
#if defined(MESPBG5P)
	  MesPbg5p(u, NumTop, rdummy, n_count);
	  cg_iter += n_count;
#endif
	}
      }
      else
	pbp_st = 0;

      /* Measure Polyakov loops */
      for(mu = 0; mu < Nd; ++mu)
	polylp (u, pollp[mu], mu);

      push(nml_out,"obsvbl1");
      Write(nml_out, TotalTrj);
      Write(nml_out, w_plaq);
      Write(nml_out, s_plaq);
      Write(nml_out, t_plaq);
      Write(nml_out, link);
      Write(nml_out, pollp);
      pop(nml_out);
      push(nml_out,"obsvbl2");
      Write(nml_out, pbp_st);
      Write(nml_out, cg_iter);
      Write(nml_out, DelH);
      Write(nml_out, DelKe);
      Write(nml_out, DelPe);
      Write(nml_out, DelFe);
      pop(nml_out);
      push(nml_out,"obsvbl3");
      Write(nml_out, ke);
      Write(nml_out, pe);
      Write(nml_out, fe);
      pop(nml_out);

      /*      if ( ( PolyEvolP == YES || RatEvolP == YES ) 
	      && FermTypeP == STAGGERED_FERMIONS && FermiP==YES )*/
      if ( FermTypeP == STAGGERED_FERMIONS && FermiP == YES ) 
	MesEV (u, KappaMC, lam_lo, lam_hi, 1);

      if ( GlueP == YES )
	fuzglue (u, j_decay, BlkAccu, BlkMax);
      if ( WilsonP == YES )
	wilslp (u, j_decay, 3);

      /* Print statistics every PrtMes measurements */
      if ( (TotalTrj % (MesTrj*PrtMes)) == 0 )
      {
	if (AlgLPStp == YES || MonitordH == YES)
          if (R0algP == NO && Algorithm != PHMCN)
	  {
	    FPRINTF(trm_out,"     <min(1,exp(-dH))> =%14.7g     <acceptance> =%14.7g\n",
		    SumAccProb/TO_REAL(n_loop),TO_REAL(n_acc)/TO_REAL(n_loop) );
	  } 
          else
	  {           
	    FPRINTF(trm_out,"     <accept-prob> =%14.7g     <acceptance> =%14.7g\n",     SumAccProb/TO_REAL(n_loop),TO_REAL(n_acc)/TO_REAL(n_loop) );

	    FPRINTF(trm_out,"     <N_Terms_Exp> =%14.7g     <N_Terms_Log> =%14.7g\n",TO_REAL(NTermsExp)/TO_REAL(n_loop),TO_REAL(NTermsLog)/TO_REAL(NTermsExp) );
	  }
      } /* end of on-line statistics display */
    } /* end of measurements */
  
    /* Calculate topological charge and the action/continuum instanton action ratio. */
    /* Calculate meson and baryon propagators only if the configuration was accepted. */
    /* For now, we do not measure psi-bar-psi from the point source since it gives */
    /* so little information. */
    if ( (TotalTrj % SpecTrj) == 0 )
    {
      /* Smeared Wilson loops (for heavy-quark potential) */
      if ( SmwlsP == YES )
      {
		
	u_smr = u;
	for(cb = 0; cb < Nsubl; ++cb)
	  u_tmp[j_decay][cb] = u_smr[j_decay][cb];

	/* Smear the space-like links numb_sm times */
	for(i = 0; i < numb_sm; ++i)
	{
	  for(mu = 0; mu < Nd; ++mu)
	    if (mu != j_decay)
	      for(cb = 0; cb < Nsubl; ++cb)
		smear (u_smr, u_tmp[mu][cb], cb, mu, 0, fact_sm, BlkAccu, BlkMax, j_decay);

	  u_smr = u_tmp;
	}

	
	push(nml_out,"Smeared_Wilson_Loops");
	Write(nml_out, numb_sm);
	Write(nml_out, fact_sm);
	pop(nml_out);

	wilslp (u_smr, j_decay, 6);
      }

      /* Topology */
      if( TopolP == YES  )
      {
	u_tmp = u;
	topol (u_tmp, TopAccu, ActAccu, NumCool, NumTop);
      }

      /* Hadron spectroscopy */
      ncg_had = 0;
      nrl_gf  = 0;

      if( MesonP == YES  )
      {
	if( spec_acc == YES )
	{
	  push(nml_out,"Do_hadron_measurements");
	  Write(nml_out, spec_acc);
	  pop(nml_out);
  
	  u_tmp = u;

	  t_srce_tmp = t_srce;
  
	  length = nrow[j_decay];
	  half_length = length/2;
  
	  t0 = t_srce[j_decay];
	  itmp1 = t0 + half_length;
	  t_srce_tmp[j_decay] = mod(itmp1, length);
  
	  /* Gauge-fix to the coulomb gauge */
	  gfix (u_tmp, j_decay, GFAccu, GFMax, nrl_gf, OrDo, OrPara);
  
	  /* Turn boundary phases on, so that they are turned off in hadron */
	  bc_spec = 1;
	  switch (FermTypeP)
	  {
	  case WILSON_FERMIONS:
	    phfctr (u_tmp, FORWARD);
	    break;
	  case STAGGERED_FERMIONS:
	    u_tmp = u_tmp * bc_phases;
	    break;
	  default:
	    QDP_error_exit("Unknown fermion type", FermTypeP);
	  }
  
	  /* For statistics, call hadron twice */
	  numKappa = 1;
	  Kappa[0] = KappaMC;
	  hadron (u_tmp, t_srce, j_decay, Kappa, numKappa, RsdCGMC, OPTION[WALL_SOURCE], OPTION[POINT_AND_WALL_SINK], ncg_had, bc_spec);
	  hadron (u_tmp, t_srce_tmp, j_decay, Kappa, numKappa, RsdCGMC, OPTION[WALL_SOURCE], OPTION[POINT_AND_WALL_SINK], ncg_had, bc_spec);
  
	  	  	    
	  push(nml_out,"cg_and_gf_iterations");
	  Write(nml_out, nrl_gf);
	  Write(nml_out, ncg_had);
	  pop(nml_out);
	  MesonP = 1;		/* spectroscopy will now only be done on accep. */
	} /* end of spectroscopy */
	else
	{
	  push(nml_out,"Do_hadron_measurements");
	  Write(nml_out, spec_acc);
	  pop(nml_out);
	  push(nml_out,"cg_and_gf_iterations");
	  Write(nml_out, nrl_gf);
	  Write(nml_out, ncg_had);
	  pop(nml_out);
	}
      }

      /* Schroedinger PCAC measurement (for FermTypeP==WILSON_FERMIONS) */
      if ( SchrPcacP == YES && SchrFun > 0 && FermTypeP == WILSON_FERMIONS )
      {
        ncg_had = 0;
	numKappa = 1;
	Kappa[0] = KappaMC;
	ClCoeff[0] = ClovCoeff;
        SFpcac (u, j_decay, Kappa, ClCoeff, ClCoeff, numKappa, RsdCGMC, ncg_had, NO, NO, 0, 0);
        push(nml_out,"PCAC_Inversion_steps");
	Write(nml_out, ncg_had);
	pop(nml_out); 
      }

      spec_acc = NO;		/* Regardless, have attempted spectroscopy */
    } /* end of spectroscopy measurements */
  
    /* This is absolutely the last namelist group Printed in this trajectory */
    push(nml_out,"EndTrajectory");
    Write(nml_out, TrjNumber);
    pop(nml_out);

    NumCG = NumCG + cg_iter;
    conf_acc = NO;		/* Be careful about placement of this line. */
  } /* end loop over trajectories */

          
  push(nml_out,"Acceptances");
  Write(nml_out, n_acc);
  pop(nml_out);

  push(nml_out,"Probability_Estimate_Exceptions");
  Write(nml_out, NProbViLa);
  Write(nml_out, NProbViSm);
  pop(nml_out);
  
  /* Close out any other code */
  END_CODE("subroutine");;
}

  
  
