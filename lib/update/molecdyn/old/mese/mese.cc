// $Id: mese.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED"

#include "chromabase.h"


//! MESE -- measures the contributions to the Hamiltonian per bosonic d.o.f.
/*!
 * \ingroup molecdyn
 *
 */
/* u      -- gauge field (Read) */
/* p_mom  -- Conjugate Momemta (Read) */
/* chi    -- pseudo-fermion field (Read) */
/* psi    -- (M_dagM)^(-1) * chi , if PolyEvolP = NO   (Read) */
/*           not used,             if PolyEvolP = YES         */
/* Beta   -- coupling constant in the relevant Hamiltonian (MC or MD) (Read) */
/* ke     -- kinetic energy ( Write ) */
/* pe     -- potential energy ( Write ) */
/* fe     -- fermionic energy ( Write ) */
/* cg_ct  -- number of CG iterations (Write) */
/* Ncb    -- number of checkerboards  ( Read ) */
/* Npf    -- Number of pseudofermion fields (Read) */

void MesE(const multi1d<LatticeColorMatrix>& u,
	  const multi1d<LatticeColorMatrix>& p_mom,
	  const multi1d<LatticeFermion>& chi,
	  const multi1d<LatticeFermion>& psi,
	  const Real& Beta,
	  const Real& Kappa,
	  Double& ke,		        /* kinetic energy   */
	  Double& pe,		        /* potential energy */
	  Double& fe,		        /* fermionic energy */
	  int cg_ct,                    /* number of CG iterations            */
                                /* in case of rational approximations */
	  int Npf)
{
  START_CODE();
  
  cg_ct = 0;

  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  Double p_mom_sq;
  MesMom(p_mom, p_mom_sq);

  /* Accumulate the fermionic energy - prefactors will be put in later */
  MesFerm(u, chi, psi, Kappa, fe, Ncb, Npf, cg_ct);

  Double w_gimp;
#if 0
  if (GlueImp != 0)
    MesGImp (u, w_gimp);
  else
    w_gimp = 0;
#endif
    
  /* Compute kinetic and fermionic energy per gauge field d.o.f.            */
  /* P_MOM_SQ is the K.E. per link, and W_FERM(i) is the F.E. per           */
  /* pseudofermion d.o.f. (of which there are Ns*Nc per site),              */
  /* corresponding to pseudofermion field i. W_GIMP is                      */
  /* the Bosonic action (= potential energy) per degree of freedom:         */
  /* N.B., there are Nd*(Nd-1)/2 plaquettes contributing to W_PLAQ, whereas */
  /* there are Nd link variables each with Nc**2 - 1 d.o.f., thus the total */
  /* bosonic action is Nd*(Nd-1)*W_PLAQ*BetaMD/2 and the bosonic action per */
  /* d.o.f. is Nd*(Nd-1)*W_PLAQ*BetaMD)/(2*Nd*(Nc*Nc-1)) =                  */
  /* (Nd-1)*W_PLAQ*BetaMD/(2*(Nc*Nc-1)).                                    */
  
  if (Nc >= 2)
  {
    pe = (Double(Nd-1)/Double(2*(Nc*Nc-1)))*Double(Beta)*w_plaq;
    pe = pe + Double(Beta)*w_gimp/Double(Nd*(Nc*Nc-1));
    ke = p_mom_sq/Double(Nc*Nc-1);
    fe = (Double(Nc*Ns)/Double(Nd*(Nc*Nc-1))) * fe;
  }
  else if (Nc == 1)
  {
    pe = (Double(Nd-1)/Double(2))*Double(Beta)*w_plaq;
    pe = pe + Double(Beta)*w_gimp/Double(Nd);
    ke = p_mom_sq;
    fe = (Double(Ns)/Double(Nd)) * fe;
  }
  else
    QDP_error_exit("Unsupported number of colors", Nc);

  /* Put in FerFactor in R-algorithm mode */
  if (R0algP == YES || Algorithm == R)
  {
    Real FerFactor;              /* Nf/2Npf for Wilson, Nf/4Npf Staggered */
  
    if (FermTypeP == WILSON_FERMIONS)
      FerFactor = Nf/TO_REAL(2*Npf);
    else if (FermTypeP == STAGGERED_FERMIONS)
      FerFactor = Nf/TO_REAL(4*Npf);
    
    fe = Double(FerFactor) * fe;
  }

  END_CODE();
}


