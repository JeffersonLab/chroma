#include "actions/ferm/invert/quda_solvers/syssolver_cg_quda_wilson_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverCGQUDAWilsonParams::SysSolverCGQUDAWilsonParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "WilsonParams", WilsonParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);

    read(paramtop, "Delta", Delta);
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverCGQUDAWilsonParams& p)
  {
    SysSolverCGQUDAWilsonParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverCGQUDAWilsonParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "WilsonParams", p.WilsonParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "Delta", p.Delta);
    pop(xml);

  }



}
