// $Id: hybtrj.cc,v 1.1 2003-12-30 19:52:28 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE AlgETrj into params of Integ. functor"

#include "chromabase.h"

using namespace QDP;

//! Performs one standard Hybrid molecular dynamic trajectory
/*!
 * \ingroup molecdyn
 *
 */
/* u         -- gauge field ( Modify ) */
/* p_mom     -- conjugate momenta to u ( Modify ) */
/* chi       -- pseudofermion field, source for CG inversion ( Read ) */
/* psi       -- present solution of M_dagger * M psi = chi ( Modify ) */
/*              not used if PolyEvolP = NO                            */
/* n_count   -- Number of CG-iterations in one trajectory ( Write ) */
/* ke        -- initial kinetic energy ( Modify ) */
/* pe        -- initial potential energy ( Modify ) */
/* fe        -- initial fermionic energy ( Modify ) */
/* Npf       -- number of pseudofermions  ( Read ) */

void HybTrj(multi1d<LatticeColorMatrix>& u,
	    multi1d<LatticeColorMatrix>& p_mom,
	    const multi1d<LatticeFermion>& chi,
	    multi1d<LatticeFermion> psi,
	    int& n_count,
	    Double& ke,			/* Kinetic energy */
	    Double& pe,			/* Bosonic action (potential energy) */
	    Double& fe,			/* Fermionic action (fermion energy) */
	    int Npf)
{
  multi1d<LatticeFermion> old_psi(Npf);
  int nn_count;
  Real t;
  Double old_stp_ke;		/* Kinetic energy at beginning of step */
  Double old_stp_pe;		/* Bosonic action (potential energy) */
  Double old_stp_fe;		/* Fermionic action (fermion energy) */
  Real DelH;			/* Change in total energy */
  Real DelKe;			/* Change in kinetic energy */
  Real DelPe;			/* Change in potential energy */
  Real DelFe;			/* Change in fermionic energy */
  bool EndP;                /* Indicates end of the trajectory */
  int i;                   /* Scratch variable: counter */
  int ichiral;

  START_CODE("HybTrj");

  push(xml_out,"HybTrj");

  old_stp_ke = ke;
  old_stp_pe = pe;
  old_stp_fe = fe;

  old_psi = psi;

  t = 0;
  EndP = false;
  Real dtH = 0.5*dt;

  n_count = 0;

  LeapP(u, p_mom, chi, chi_norm, psi, dtH, Ncb, Npf, nn_count);
  n_count = n_count + nn_count;

  while ( ! EndP )
  {
    t += dt;

    LeapU (p_mom, u, dt);
    if (PolyEvolP && RatEvolP)      // uses global bool PolyEvolP,RatEvolP
    {
      Interpol(psi, old_psi, Real(-1), Npf);
    }
    EndP = EndOfTrj(t);
    LeapPMX(u, p_mom, chi, chi_norm, psi, dtH, dtH, nn_count, ke, pe, fe, EndP, Npf);
    n_count += nn_count;

    if ( MonitordH )   // uses global bool MonitordH
    {
      /* Find the change in energy per d.o.f. during the current leapfrog step */
      /* N.B.!! Potential energy appears with a minus sign! */
      
      DelKe = ke - old_stp_ke;
      DelPe = pe - old_stp_pe;
      DelFe = fe - old_stp_fe;
      DelH  = DelKe - DelPe + DelFe;

      push(xml_out,"HybStpMon");
      Write(xml_out, t);
      Write(xml_out, nn_count);
      Write(xml_out, DelH);
      Write(xml_out, DelKe);
      Write(xml_out, DelPe);
      Write(xml_out, DelFe);
      pop(xml_out);

      if ( FermAct == OVERLAP_POLE && FermiP && ! EndP )
	ZeroSector(u, 5, ichiral, Ncb);

      old_stp_ke = ke;
      old_stp_pe = pe;
      old_stp_fe = fe;
    }
  }

  pop(xml_out);
  
  END_CODE("HybTrj");
}

