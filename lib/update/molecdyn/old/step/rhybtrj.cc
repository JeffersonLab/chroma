// $Id: rhybtrj.cc,v 3.0 2006-04-03 04:59:11 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE GLOBAL params of Integ. functor"

#include "chromabase.h"


//! Performs one R_algorithm type trajectory (see phys.rev.D v35,no8,p2531)
/*!
 * \ingroup molecdyn
 *
 */
/* u         -- gauge field ( Modify ) */
/* p_mom     -- conjugate momenta to u ( Modify ) */
/* chi       -- pseudofermion field, source for CG inversion ( Read ) */
/* chi_norm  -- any moron would understand this ( Read ) */
/* psi       -- present solution of M_dagger * M psi = chi ( Modify ) */
/* n_count   -- Number of CG-iterations in one trajectory ( Write ) */
/* ke        -- initial kinetic energy ( Modify ) */
/* pe        -- initial potential energy ( Modify ) */
/* fe        -- initial fermionic energy ( Modify ) */
/* Npf       -- number of pseudofermions  ( Read ) */

void RHybTrj(multi1d<LatticeColorMatrix>& u,
	     multi1d<LatticeColorMatrix>& p_mom,
	     const multi1d<LatticeFermion>& chi,
	     multi1d<LatticeFermion>& psi,
	     int& n_count,
	     Double& ke,			/* Kinetic energy */
	     Double& pe,			/* Bosonic action (potential energy) */
	     Double& fe,			/* Fermionic action (fermion energy) */
	     int Npf)
{
  START_CODE();

  int nn_count;
  Real t;
  Double old_stp_ke;		/* Kinetic energy at beginning of step */
  Double old_stp_pe;		/* Bosonic action (potential energy) */
  Double old_stp_fe;		/* Fermionic action (fermion energy) */
  Real DelH;			/* Change in total energy */
  Real DelKe;			/* Change in kinetic energy */
  Real DelPe;			/* Change in potential energy */
  Real DelFe;			/* Change in fermionic energy */
  bool EndP;                    /* Indicates end of the trajectory */
  Real FerFactor;
  Real R1_dt;
  Real R2_dt;
  Real R3_dt;
  Real R4_dt;

  push(nml_out,"RHybTrj");

  old_stp_ke = ke;
  old_stp_pe = pe;
  old_stp_fe = fe;

  t = 0;
  EndP = false;

  if (FermTypeP == WILSON_FERMIONS)
    FerFactor = Nf/TO_REAL(2*Npf);
  else if (FermTypeP == STAGGERED_FERMIONS)
    FerFactor = Nf/TO_REAL(4*Npf);
  R4_dt = WORD_VALUE(WORD_dt,HALF)*dt;
  R1_dt = (1-FerFactor)*R4_dt;
  R2_dt = FerFactor*R4_dt;
  R3_dt = (TO_REAL(2)-FerFactor)*R4_dt;

  n_count = 0;

  LeapU(p_mom, u, R1_dt);
  while( ! EndP )
  {
    t += dt;
    FrmNse(u, chi, psi, chi_norm, Ncb, Ncb, Npf, nn_count);
    n_count += nn_count;

    LeapU(p_mom, u, R2_dt);

    LeapP(u, p_mom, chi, chi_norm, psi, dt, Npf, nn_count);
    n_count += nn_count;

    EndP = EndOfTrj(t);
    if (! EndP) 
      LeapU(p_mom, u, R3_dt);

    if ( MonitordH )
    {
      /* Find the change in energy per d.o.f. during the current leapfrog step*/
      /* N.B.!! Potential energy appears with a minus sign! */
      
      DelKe = ke - old_stp_ke;
      DelPe = pe - old_stp_pe;
      DelFe = fe - old_stp_fe;
      DelH  = DelKe - DelPe + DelFe;

      push(nml_out,"HybStpMon");
      write(nml_out, "t", t);
      write(nml_out, "nn_count", nn_count);
      write(nml_out, "DelH", DelH);
      write(nml_out, "DelKe", DelKe);
      write(nml_out, "DelPe", DelPe);
      write(nml_out, "DelFe", DelFe);
      pop(nml_out);

      old_stp_ke = ke;
      old_stp_pe = pe;
      old_stp_fe = fe;
    }    
  }
  LeapU(p_mom, u, R4_dt);

  pop(xml_out);

  END_CODE();
}

