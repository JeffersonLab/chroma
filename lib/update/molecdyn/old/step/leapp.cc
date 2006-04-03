// $Id: leapp.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE GLOBAL PARAMS INTO FUNCTOR"

#include "chromabase.h"


//! Leap one fractional step in the momenta
/*!
 * \ingroup molecdyn
 *
 */
/* u        -- gauge field ( Read ) */
/* p_mom    -- conjugate momenta to u ( Modify ) */
/* chi      -- pseudofermion field. ( Read ) */
/* chi_norm -- | chi | ( Read ) */
/* psi      -- Pseudofermion field. Contains the solution of M_dagger*M psi = chi */
/*             not used if PolyEvolP = YES ( Read ) */
/* eps      -- step size to use ( Read ) */
/* cg_count -- number of CG iterations in force routine ( Write ) */
/* Npf      -- number of pseudofermions  ( Read ) */

void LeapP(const multi1d<LatticeColorMatrix>& u,
	   multi1d<LatticeColorMatrix>& p_mom,
	   const multi1d<LatticeFermion>& chi,
	   const multi1d<LatticeFermion>& psi,
	   const Real& eps,
	   int Npf,
	   int& cg_count)
{
  START_CODE();

  multi1d<LatticeColorMatrix> ds_u(Nd);
  Real FerFactor;              /* Nf/2Npf for Wilson, Nf/4Npf Staggered */
  
  cg_count=0;

  dsdu(u, ds_u);
  
  Real dummy = BetaMD/Real(2*Nc);
  if (R0algP || Algorithm == R)    // global params   bool R0algP, int Algorithm
  {
    if (FermTypeP == WILSON_FERMIONS)
      FerFactor = Nf/Real(2*Npf);
    else if (FermTypeP == STAGGERED_FERMIONS)
      FerFactor = Nf/Real(4*Npf);

    dummy /= FerFactor;
  }
  /* Note that the 1/FerFactor will be removed from this gauge    */
  /* part contribution once the fermionic part is included.       */

  for(int mu=0; mu<Nd; ++mu)
    ds_u[mu] *= dummy;
  
  if ( FermiP )    // global bool FermiP
  {
    phfctr(u, FORWARD);		        /* ON */
    if (PolyEvolP)            // global bool PolyEvolP, RatEvolP
    {
      for(int i = 0; i < Npf; ++i)
	polydsduf(u, ds_u, chi[i][0]);
    }
    else if (RatEvolP == YES)
    {
      LatticeFermion tmp;

      /* tmp could be used for energy measurement, but  */
      /* probably only in pqp leapfrog scheme. Here not */
      /* used.                                          */     
      for(int i = 0; i < Npf; ++i)
      {
	int cg_ct;

	ratdsduf(u, ds_u, chi[i][0], tmp, cg_ct);
	cg_count += cg_ct;
      }
    }
    else
    {
      dsduf(u, ds_u, chi, chi_norm, psi, Npf, cg_count);
    }
    phfctr(u, BACKWARD);		        /* OFF */
  }

#if 0
  /* If using Schroedinger functional, zero out the boundaries */
  if ( SchrFun > 0 )
  {
    FILLMASK(ds_u, lSFmask, ZERO);
  }
#endif

  /* Put the FerFactor in if in R0 or R modes */
  if (R0algP == YES || Algorithm == R)
  {
    for(int mu=0; mu<Nd; ++mu)
      ds_u[mu] *= FerFactor;
  }

  /* p_mom = p_mom - eps*[Ds/Du + DS_f/Du] */
  for(int mu=0; mu<Nd; ++mu)
    p_mom[mu] -= ds_u[mu] * eps;
  
    
  for(int mu=0; mu<Nd; ++mu)
    taproj(p_mom[mu]);
  
  END_CODE();
}
