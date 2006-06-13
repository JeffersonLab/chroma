// $Id: szin_io.cc,v 3.1 2006-06-13 18:18:14 bjoo Exp $

/*! \file
 *  \brief Reader/writers for szin headers
 */

#include "chromabase.h"
#include "chroma_config.h"
#include "io/szin_io.h"
#include <time.h>

namespace Chroma 
{

  // Initialize header with default values
  SzinGauge_t::SzinGauge_t()
  {
    cfg_version = 7;

    FermTypeP = 0;
    Nd = QDP::Nd;
    Nc = QDP::Nc;
    BetaMC = 0;
    BetaMD = 0;

    KappaMC = 0;
    KappaMD = 0;
    MassMC = 0;
    MassMD = 0;
    dt = 0;
    MesTrj = 0;
    TotalCG = 0;
    TotalTrj = 0;
    spec_acc = 1;

    NOver = 0;
    TotalTry = 0;
    TotalFail = 0;
    Nf = 0;
    Npf = 0;
    RefMomTrj = 0;
    RefFnoiseTrj = 0;
    LamPl = 0;
    LamMi = 0;
    AlpLog = 0;
    AlpExp = 0;

    nrow = Layout::lattSize();
    RNG::savern(seed);

    banner = CHROMA_PACKAGE_STRING;
    banner += ", ";
    banner += QDP_PACKAGE_STRING;

    {
      time_t now;

      if(time(&now)==-1)
      {
	QDPIO::cerr<<"SzinGauge_t: cannot get the time.\n";
	QDP_abort(1);
      }
      tm *tp = localtime(&now);

      char date_tmp[64];
      strftime(date_tmp, 63, "%d %b %y %X %Z", tp);
      date = date_tmp;
    }
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

}  // end namespace Chroma
