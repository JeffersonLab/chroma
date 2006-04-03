// $Id: leappmx.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE GLOBAL PARAMS INTO FUNCTOR"

#include "chromabase.h"


//! Leaps P forward step eps with optional monitor of DelH
/*!
 * \ingroup molecdyn
 *
 */
/* The trajectory is always carried out with CG residual RsdCGMD.    */
/* If MonitordH is set and RsdCGMDMax > RsdCGMD and PolyEvolP = NO   */
/* then the change in DelH is monitored as the residual is decreased */
/* by factor(s) of RsdCGMDFctr from RsdCGMDMax to RsdCGMD. */

/* u         -- gauge field ( Modify ) */
/* p_mom     -- conjugate momenta to u ( Modify ) */
/* chi       -- pseudofermion field, source for CG inversion ( Read ) */
/* chi_norm  -- any moron would understand this ( Read ) */
/* psi       -- present solution of M_dagger * M psi = chi ( Modify ) */
/*              not modified if PolyEvolP = YES   */ 
/* eps1      -- first half-step size ( Read ) */
/* eps2      -- second half-step size ( Read ) */
/* n_count   -- Number of CG-iterations in one trajectory ( Write ) */
/* ke        -- kinetic energy ( Write ) */
/* pe        -- potential energy ( Write ) */
/* fe        -- fermionic energy ( Write ) */
/* EndP      -- Flag for exit after P half step ( Read ) */
/* Npf       -- number of pseudofermions  ( Read ) */

void LeapPMX(multi1d<LatticeColorMatrix>& u,
	     const multi1d<LatticeColorMatrix>& p_mom,
	     const multi1d<LatticeFermion>& chi,
	     const multi1d<LatticeFermion>& psi,
	     const Real& eps1,
	     const Real& eps2,
	     bool EndP,
	     int& n_count,
	     Double& ke,		        /* kinetic energy */
	     Double& pe,		        /* potential energy */
	     Double& fe,		        /* fermionic energy */
	     int Npf)
{
  START_CODE();

  Real RsdCG;                  /* Temporary for residual */
  Double ke_t;		      /* kinetic energy for least accurate inv. */
  Double pe_t;		      /* potential energy for least accurate inversion */
  Double fe_t;		      /* fermionic energy for least accurate inversion */
  Real DelHCG;	              /* Change in total energy for least accurate inversion */
  Real DelKeCG;	      /* Change in kinetic energy for least accurate inversion */
  /*Real DelPeCG;		# Change in potential energy for least accurate inversion */
  Real DelFeCG;		/* Change in fermionic energy for least accurate inversion */
  int del_n_count;         /* CG iterations for each incremental step */
  int nn_count;            /* CG iterations in leapp or mese         */
  Real eps12;                  /* Scratch variable */
  multi1d<Double> psi_norm(Npf);   // only for monitor output

  push(xml_out,"LeapPMX");

  eps12 = eps1 + eps2;

  /* Optimization */
  n_count = 0;
  if ( MonitordH == NO && AlgLPStp == NO && EndP == YES && RefNextTrj == YES ) 
  {    // lots of global params here
    END_CODE();
    return;
  }

  if ( FermiP && ! PolyEvolP && ! RatEvolP )  // more global params
  {
    RsdCG = max(RsdCGMD, RsdCGMDMax);

    if ( ( FermAct == OVERLAP_POLE || FermAct == TRUNC_OVERLAP )
	 && NOperEigDim > 0 )
      CalcEigVec (u, NO);

    phfctr (u, FORWARD);              /* ON */
    for(int i = 0; i < Npf; ++i)
    {
      invert(u, chi[i], psi[i], KappaMD, RsdCG, del_n_count);
      n_count += del_n_count;
    }
    phfctr (u, BACKWARD);             /* OFF */

    if ( MonitordH == YES )
    {
      for(int i = 0; i < Npf; ++i)
	psi_norm[i] = sqrt(norm2(psi[i]));
    }
  }


  if ( ! MonitordH && ! EndP )
  {
    LeapP(u, p_mom, chi, psi, eps12, Npf, nn_count);
    n_count += nn_count;
  }
  else
  {
    if ( RsdCGMD < RsdCGMDMax && PolyEvolP == NO && RatEvolP == NO)
    {
      multi1d<LatticeColorMatrix> p_mom_t(Nd);  /* temp. mom. for inc. steps */

      /* Calculate the energies at the least accurate residual */
      MesE (u, p_mom, chi, psi, BetaMD, KappaMD, ke_t, pe_t, fe_t, Npf, nn_count);
      n_count = n_count + nn_count;

      while(1)
      {
	p_mom_t = p_mom;
	LeapP (u, p_mom_t, chi, psi, eps1, Npf, nn_count);
        n_count = n_count + nn_count;
         
	/* Calculate the energies for more accurate inversions */
	MesE(u, p_mom_t, chi, psi, BetaMD, KappaMD, ke, pe, fe, Npf, nn_count);
        n_count = n_count + nn_count;

	/* Find the change in energy per d.o.f. during the current leapfrog step */
	/* N.B.!! Potential energy appears with a minus sign! */
	/*     DelPeCG = pe - pe_t = 0; */
	DelKeCG = ke - ke_t;
	DelFeCG = fe - fe_t;
	DelHCG  = DelKeCG + DelFeCG;

	push(xml_out,"RsdCGMon");
	write(xml_out, "del_n_count", del_n_count);
	write(xml_out, "n_count", n_count);
	write(xml_out, "psi_norm", psi_norm);
	write(xml_out, "RsdCG",  RsdCG);
	write(xml_out, "DelHCG", DelHCG);
	write(xml_out, "DelKeCG", DelKeCG);
	write(xml_out, "DelFeCG", DelFeCG);
	pop(xml_out);

	if ( RsdCG <= RsdCGMD )
	{
	  p_mom = p_mom_t;
	  break;
	}
	else
	{
	  RsdCG = max(RsdCG*RsdCGMDFctr, RsdCGMD);
	  
	  phfctr(u, FORWARD);              /* ON */
          for( i = 0; i < Npf; ++i )
          {
	    invert(u, chi[i], psi[i], KappaMD, RsdCG, del_n_count);
	    n_count += del_n_count;
          }
	  phfctr(u, BACKWARD);             /* OFF */

          for( i = 0; i < Npf; ++i )
          {
	    psi_norm[i] = sqrt(norm2(psi[i]));
          }
	} /* end if RsdCG */
      } /* end while */
    } /* end then  */
    else if (RsdCGMD >= RsdCGMDMax || PolyEvolP == YES || RatEvolP == YES)
    {
      LeapP(u, p_mom, chi, psi, eps1, Npf, nn_count);
      n_count = n_count + nn_count;

      if ( MonitordH == YES ) 
      {
/* Calculate the Hamiltonian at the end of a leapfrog step */
	MesE(u, p_mom, chi, psi, BetaMD, KappaMD, ke, pe, fe, Npf, nn_count);
        n_count += nn_count;

      }
    } /* end else */

    if ( EndP ) 
    {
      END_CODE();
      return;
    }

    LeapP(u, p_mom, chi, psi, eps2, Npf, nn_count);
    n_count += nn_count;
  }

  pop(xml_out);

  END_CODE();
}
