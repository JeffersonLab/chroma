// $Id: szin_io.cc,v 1.4 2004-02-23 03:08:02 edwards Exp $

/*! \file
 *  \brief Reader/writers for szin headers
 */

#include "chromabase.h"
#include "io/szin_io.h"

using namespace QDP;

// Initialize header with default values
void initHeader(SzinGauge_t& header)
{
  header.cfg_version = 7;

  header.FermTypeP = 0;
  header.Nd = Nd;
  header.Nc = Nc;
  header.BetaMC = 0;
  header.BetaMD = 0;

  header.KappaMC = 0;
  header.KappaMD = 0;
  header.MassMC = 0;
  header.MassMD = 0;
  header.dt = 0;
  header.MesTrj = 0;
  header.TotalCG = 0;
  header.TotalTrj = 0;
  header.spec_acc = 1;

  header.NOver = 0;
  header.TotalTry = 0;
  header.TotalFail = 0;
  header.Nf = 0;
  header.Npf = 0;
  header.RefMomTrj = 0;
  header.RefFnoiseTrj = 0;
  header.LamPl = 0;
  header.LamMi = 0;
  header.AlpLog = 0;
  header.AlpExp = 0;

  header.nrow = Layout::lattSize();
  RNG::savern(header.seed);

  header.banner = "chroma stuff - need more here";
  header.date = "put in some date here";
}



// Source header read
void read(XMLReader& xml, const string& path, SzinGauge_t& header)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "cfg_version", header.cfg_version);

  read(paramtop, "date", header.date);
  read(paramtop, "banner", header.banner);

  read(paramtop, "FermTypeP", header.FermTypeP);
  read(paramtop, "Nd", header.Nd);
  read(paramtop, "Nc", header.Nc);
  read(paramtop, "BetaMC", header.BetaMC);
  read(paramtop, "BetaMD", header.BetaMD);

  read(paramtop, "KappaMC", header.KappaMC);
  read(paramtop, "KappaMD", header.KappaMD);
  read(paramtop, "MassMC", header.MassMC);
  read(paramtop, "MassMD", header.MassMD);
  read(paramtop, "dt", header.dt);
  read(paramtop, "MesTrj", header.MesTrj);
  read(paramtop, "TotalCG", header.TotalCG);
  read(paramtop, "TotalTrj", header.TotalTrj);
  read(paramtop, "spec_acc", header.spec_acc);

  read(paramtop, "NOver", header.NOver);
  read(paramtop, "TotalTry", header.TotalTry);
  read(paramtop, "TotalFail", header.TotalFail);
  read(paramtop, "Nf", header.Nf);
  read(paramtop, "Npf", header.Npf);
  read(paramtop, "RefMomTrj", header.RefMomTrj);
  read(paramtop, "RefFnoiseTrj", header.RefFnoiseTrj);
  read(paramtop, "LamPl", header.LamPl);
  read(paramtop, "LamMi", header.LamMi);
  read(paramtop, "AlpLog", header.AlpLog);
  read(paramtop, "AlpExp", header.AlpExp);

  read(paramtop, "nrow", header.nrow);
  read(paramtop, "seed", header.seed);
}


//! Source header writer
void write(XMLWriter& xml, const string& path, const SzinGauge_t& header)
{
  push(xml, path);

  write(xml, "cfg_version", header.cfg_version);

  write(xml, "date", header.date);
  write(xml, "banner", header.banner);

  write(xml, "FermTypeP", header.FermTypeP);
  write(xml, "Nd", header.Nd);
  write(xml, "Nc", header.Nc);
  write(xml, "BetaMC", header.BetaMC);
  write(xml, "BetaMD", header.BetaMD);

  write(xml, "KappaMC", header.KappaMC);
  write(xml, "KappaMD", header.KappaMD);
  write(xml, "MassMC", header.MassMC);
  write(xml, "MassMD", header.MassMD);
  write(xml, "dt", header.dt);
  write(xml, "MesTrj", header.MesTrj);
  write(xml, "TotalCG", header.TotalCG);
  write(xml, "TotalTrj", header.TotalTrj);
  write(xml, "spec_acc", header.spec_acc);

  write(xml, "NOver", header.NOver);
  write(xml, "TotalTry", header.TotalTry);
  write(xml, "TotalFail", header.TotalFail);
  write(xml, "Nf", header.Nf);
  write(xml, "Npf", header.Npf);
  write(xml, "RefMomTrj", header.RefMomTrj);
  write(xml, "RefFnoiseTrj", header.RefFnoiseTrj);
  write(xml, "LamPl", header.LamPl);
  write(xml, "LamMi", header.LamMi);
  write(xml, "AlpLog", header.AlpLog);
  write(xml, "AlpExp", header.AlpExp);

  write(xml, "nrow", header.nrow);
  write(xml, "seed", header.seed);

  pop(xml);
}

