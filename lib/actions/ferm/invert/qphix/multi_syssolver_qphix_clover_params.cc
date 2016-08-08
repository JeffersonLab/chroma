#include "actions/ferm/invert/qphix/multi_syssolver_qphix_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {


  MultiSysSolverQPhiXCloverParams::MultiSysSolverQPhiXCloverParams(XMLReader& xml, 
							 const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "MaxShifts", MaxShifts);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", CloverParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);
    
    if ( paramtop.count("Verbose") > 0 ) { 
      read(paramtop, "Verbose", VerboseP);
    }
    else { 
      VerboseP = false;
    }

    if ( paramtop.count("SolutionCheck") > 0 ) { 
      read(paramtop, "SolutionCheck", SolutionCheckP);
    }
    else { 
      SolutionCheckP = true;
    }

    if (paramtop.count("Delta") > 0 ) { 
      read(paramtop,"Delta", Delta);
    }
    else { 
      Delta = Real(-1);
    }


    if( paramtop.count("RsdToleranceFactor") > 0 ) { 
      read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
    }
    else { 
      RsdToleranceFactor = Real(10); // Tolerate an order of magnitude difference by default.
    }

    if( paramtop.count("Tune") > 0 ) { 
      read(paramtop, "Tune", TuneP);
    }
    else { 
      TuneP = false;
    }
    
  }
  
  void read(XMLReader& xml, const std::string& path, 
	    MultiSysSolverQPhiXCloverParams& p)
  {
    MultiSysSolverQPhiXCloverParams tmp(xml, path);
    p = tmp;
  }
  
  void write(XMLWriter& xml, const std::string& path, 
	     const MultiSysSolverQPhiXCloverParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);

    // Hack. Write delta only if it is greater than 0. 
    // -ve sign indicates user did not supply
    if( toBool(p.Delta > 0 ) ) {
      write(xml, "Delta", p.Delta);
    }
    write(xml, "MaxShifts", p.MaxShifts);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "Verbose", p.VerboseP);
    write(xml, "SolutionCheck", p.SolutionCheckP);
    write(xml, "CloverParams", p.CloverParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);
    write(xml, "Tune", p.TuneP);
    

  }



}
