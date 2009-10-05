#include "actions/ferm/invert/quda_solvers/syssolver_quda_wilson_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverQUDAWilsonParams::SysSolverQUDAWilsonParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "WilsonParams", WilsonParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);

    read(paramtop, "Delta", Delta);
    read(paramtop, "SolverType", solverType);
    if (solverType != "CG" && solverType != "BICGSTAB" ) { 
      QDPIO::cout << "Supported solver types are CG and BICGSTAB" << endl;
      QDPIO::cout << "You entered " << solverType << endl;
      QDP_abort(1);
    }
    if ( paramtop.count("Verbose") > 0 ) { 
      read(paramtop, "Verbose", verboseP);
    }
    else { 
      verboseP = false;
    }

  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverQUDAWilsonParams& p)
  {
    SysSolverQUDAWilsonParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverQUDAWilsonParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "WilsonParams", p.WilsonParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "Delta", p.Delta);
    write(xml, "SolverType", p.solverType);
    write(xml, "Verbose", p.verboseP);
      

    pop(xml);

  }



}
