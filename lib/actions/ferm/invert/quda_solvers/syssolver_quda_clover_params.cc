#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverQUDACloverParams::SysSolverQUDACloverParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", CloverParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);

    read(paramtop, "Delta", Delta);
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverQUDACloverParams& p)
  {
    SysSolverQUDACloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverQUDACloverParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.CloverParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "Delta", p.Delta);
    pop(xml);

  }



}
