// $Id: szin_io.h,v 1.3 2004-02-13 14:59:04 bjoo Exp $

/*! \file
 *  \brief  Routines associated with SZIN gauge field IO
 */

#ifndef __szin_io_h__
#define __szin_io_h__

#include <string>

//! Szin gauge field header
struct SzinGauge_t
{
  multi1d<int> nrow;    // Lattice size
  int Nd;               // Number of spacetime dimensions
  int Nc;               // Number of colors

  int TotalTrj;         // Total number of trajectories
  int TotalCG;          // Total number of CG iterations
  int FermTypeP;        // Fermion type
  int spec_acc;         /* Acceptance flag for spectroscopy */
//int MesItr;		/* Iterations per measurement */
//int TotalItr;		/* Total number of iterations */
  int NOver;		/* Number of overrelaxation steps */
  int TotalTry;		/* Total number of heatbath trials */
  int TotalFail;	/* Total number of heatbath failures */
  int Npf;              /* Number of pseudofermions */
  int RefMomTrj;        /* Trajectories per momentum refreshment */
  int RefFnoiseTrj;     /* Trajectories per fermion-noise refreshment */
  Real32 MesTrj;        // Trajectories per measurement (as a Float)
  Real32 BetaMC;	/* 6/g**2 */
  Real32 BetaMD;	/* 6/g**2 (MD) */
  Real32 dt;		/* Step size */
  Real32 KappaMC;	/* Hopping parameter */
  Real32 KappaMD;	/* Hopping parameter */
  Real32 MassMC;        /* MC mass */
  Real32 MassMD;        /* MC mass */
  Real32 Nf;            /* Number of flavours */
  Real32 LamPl;         /* Stochastic acc/rej parameter */
  Real32 LamMi;         /* Stochastic acc/rej parameter */
  Real32 AlpLog;        /* For estimate of Log(1+x)     */
  Real32 AlpExp;        /* For estimate of Exp(x)       */
  QDP::Seed   seed;		/* Random number seed */

  int cfg_version;      /* Configuration file version number */

  std::string banner;
  std::string date;
};


//! Initialize header with default values
void szinGaugeInit(SzinGauge_t& header);

//! Source header read
void read(XMLReader& xml, const std::string& path, SzinGauge_t& header);

//! Source header writer
void write(XMLWriter& xml, const std::string& path, const SzinGauge_t& header);

#endif
