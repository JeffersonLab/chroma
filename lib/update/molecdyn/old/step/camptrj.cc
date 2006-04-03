// $Id: camptrj.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

#error "NOT FULLY CONVERTED - NEED TO MOVE AlgETrj into params of Integ. functor"

#include "chromabase.h"


//! Performs one standard Hybrid molecular dynamic trajectory.
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

void campTrj(multi1d<LatticeColorMatrix> u,
	     multi1d<LatticeColorMatrix> p_mom,
	     multi1d<LatticeFermion> chi(Npf),
	     multi1d<LatticeFermion> psi(Npf),
	     int n_count,
	     multi1d<Double> chi_norm(Npf),
	     Double& ke,			/* Kinetic energy */
	     Double& pe,			/* Bosonic action (potential energy) */
	     Double& fe,			/* Fermionic action (fermion energy) */
	     int Npf)
{
  START_CODE();
  
  multi1d<LatticeFermion> old_psi(Npf);
  Real sigma;			/* This is 2**(1/3) */
  int nn_count;
  Real t;
  Double old_stp_ke;		/* Kinetic energy at beginning of step */
  Double old_stp_pe;		/* Bosonic action (potential energy) */
  Double old_stp_fe;		/* Fermionic action (fermion energy) */
  Double old_wgl_ke;		/* Kinetic energy at beginning of wiggle */
  Double old_wgl_pe;		/* Bosonic action (potential energy) */
  Double old_wgl_fe;		/* Fermionic action (fermion energy) */
  Real DelH;			/* Change in total energy */
  Real DelKe;			/* Change in kinetic energy */
  Real DelPe;			/* Change in potential energy */
  Real DelFe;			/* Change in fermionic energy */
  Real eps;
  Real epsH;
  Real epsUBac;
  Real epsSigH;
  Real sigmainv;
  Real minus1;
  int wgl_n_count;
  bool EndP                 /* Flag indicating end of the trajectory */
  int i;                    /* Scratch variable: counter */
  
  push(xml_out,"campTrj");

  old_wgl_ke = ke;
  old_wgl_pe = pe;
  old_wgl_fe = fe;
  
  old_stp_ke = ke;
  old_stp_pe = pe;
  old_stp_fe = fe;
  
  old_psi = psi;
  
  t = 0;
  EndP = false;
  
  n_count = 0;
  
  /* Step sizes for a campostrini step */
  sigma   =  pow(TO_REAL(2),(TO_REAL(1)/TO_REAL(3)));
  eps     =  dt / (TO_REAL(2) - sigma);
  epsH    =  WORD_VALUE(WORD_eps,HALF) * eps;
  epsUBac = -sigma * eps;
  epsSigH = -sigma * epsH;
  sigmainv = TO_REAL(1) / sigma;
  minus1  = TO_REAL(-1);
  
  LeapP(u, p_mom, chi, chi_norm, psi, epsH, Npf, nn_count);
  
  n_count += nn_count;

  while ( ! EndP )
  {
    t = t + dt;
    wgl_n_count = 0;

    /* First step of wiggle */
    LeapU(p_mom, u, eps);

    Interpol(psi, old_psi, minus1, Npf);

    LeapPMX (u, p_mom, chi, chi_norm, psi, epsH, epsSigH, nn_count, ke, pe, fe, NO, Npf);
    n_count = n_count + nn_count;

    if ( MonitordH )
    {
      wgl_n_count = wgl_n_count + nn_count;


      /* Find the change in energy per d.o.f. during the current leapfrog step */
      /* N.B.!! Potential energy appears with a minus sign! */

      DelKe = ke - old_stp_ke;
      DelPe = pe - old_stp_pe;
      DelFe = fe - old_stp_fe;
      DelH  = DelKe - DelPe + DelFe;

      push(xml_out,"CampStpMon1");
      write(xml_out, "t", t);
      write(xml_out, "nn_count", nn_count);
      write(xml_out, "DelH", DelH);
      write(xml_out, "DelKe", DelKe);
      write(xml_out, "DelPe", DelPe);
      write(xml_out, "DelFe", DelFe);
      pop(xml_out);

      old_stp_ke = ke;
      old_stp_pe = pe;
      old_stp_fe = fe;
    }

    /* Second step of wiggle */
    /*_ */
    LeapU(p_mom, u, epsUBac);

    Interpol(psi, old_psi, sigma, Npf);

    LeapPMX(u, p_mom, chi, chi_norm, psi, epsSigH, epsH, nn_count, ke, pe, fe, NO, Npf);
    n_count = n_count + nn_count;

    if ( MonitordH )
    {
      wgl_n_count = wgl_n_count + nn_count;


      /* Find the change in energy per d.o.f. during the current leapfrog step */
      /* N.B.!! Potential energy appears with a minus sign! */

      DelKe = ke - old_stp_ke;
      DelPe = pe - old_stp_pe;
      DelFe = fe - old_stp_fe;
      DelH  = DelKe - DelPe + DelFe;

      push(xml_out,"CampStpMon2");
      write(xml_out, "t", t);
      write(xml_out, "nn_count", nn_count);
      write(xml_out, "DelH", DelH);
      write(xml_out, "DelKe", DelKe);
      write(xml_out, "DelPe", DelPe);
      write(xml_out, "DelFe", DelFe);
      pop(xml_out);

      old_stp_ke = ke;
      old_stp_pe = pe;
      old_stp_fe = fe;
    }

    /* Third step of wiggle */
    LeapU(p_mom, u, eps);

    Interpol(psi, old_psi, sigmainv, Npf);

    EndP = EndOfTrj(t);

    LeapPMX(u, p_mom, chi, chi_norm, psi, epsH, epsH, nn_count, ke, pe, fe, EndP, Npf);
    n_count = n_count + nn_count;

    if ( MonitordH )
    {
      wgl_n_count = wgl_n_count + nn_count;


      /* Find the change in energy per d.o.f. during the current leapfrog step */
      /* N.B.!! Potential energy appears with a minus sign! */

      DelKe = ke - old_stp_ke;
      DelPe = pe - old_stp_pe;
      DelFe = fe - old_stp_fe;
      DelH  = DelKe - DelPe + DelFe;

      push(xml_out,"CampStpMon3");
      write(xml_out, "t", t);
      write(xml_out, "nn_count", nn_count);
      write(xml_out, "DelH", DelH);
      write(xml_out, "DelKe", DelKe);
      write(xml_out, "DelPe", DelPe);
      write(xml_out, "DelFe", DelFe);
      pop(xml_out);

      old_stp_ke = ke;
      old_stp_pe = pe;
      old_stp_fe = fe;


      /* Find the change in energy per d.o.f. during the current Campostrini wiggle */

      DelKe = ke - old_wgl_ke;
      DelPe = pe - old_wgl_pe;
      DelFe = fe - old_wgl_fe;
      DelH  = DelKe - DelPe + DelFe;

      push(xml_out,"CampWglMon");
      write(xml_out, "t", t);
      write(xml_out, "wgl_n_count", wgl_n_count);
      write(xml_out, "DelH", DelH);
      write(xml_out, "DelKe", DelKe);
      write(xml_out, "DelPe", DelPe);
      write(xml_out, "DelFe", DelFe);
      pop(xml_out);

      old_wgl_ke = ke;
      old_wgl_pe = pe;
      old_wgl_fe = fe;
    }
  }
  
  pop(xml_out);

  END_CODE();
}
