#include "actions/ferm/invert/syssolver_cg_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverCGCloverParams::SysSolverCGCloverParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "MaxCG", MaxCG);
    read(paramtop, "RsdCG", RsdCG);
    read(paramtop, "CloverParams", clovParams);
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverCGCloverParams& p)
  {
    SysSolverCGCloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverCGCloverParams& p) {
    push(xml, path);
    write(xml, "MaxCG", p.MaxCG);
    write(xml, "RsdCG", p.RsdCG);
    write(xml, "CloverParams", p.clovParams);
    pop(xml);

  }



}
