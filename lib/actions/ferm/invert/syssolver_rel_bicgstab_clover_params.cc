#include "actions/ferm/invert/syssolver_rel_bicgstab_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverReliableBiCGStabCloverParams::SysSolverReliableBiCGStabCloverParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", clovParams);
    read(paramtop, "Delta", Delta);
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverReliableBiCGStabCloverParams& p)
  {
    SysSolverReliableBiCGStabCloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverReliableBiCGStabCloverParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.clovParams);
    write(xml, "Delta", p.Delta);
    pop(xml);

  }



}
