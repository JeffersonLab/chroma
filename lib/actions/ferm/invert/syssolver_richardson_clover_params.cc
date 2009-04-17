#include "actions/ferm/invert/syssolver_richardson_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverRichardsonCloverParams::SysSolverRichardsonCloverParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", clovParams);
    innerSolverParams = readXMLGroup(paramtop, "InnerSolverParams", "invType");
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverRichardsonCloverParams& p)
  {
    SysSolverRichardsonCloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverRichardsonCloverParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.clovParams);
    xml << p.innerSolverParams.xml;
    pop(xml);

  }



}
