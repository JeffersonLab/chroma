#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/quda_gcr_params.h"

using namespace QDP;

namespace Chroma {
  GCRInnerSolverParams::GCRInnerSolverParams(XMLReader& xml, 
					     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "RsdSloppy", tolSloppy);
    read(paramtop, "MaxIterSloppy", maxIterSloppy);
    read(paramtop, "NKrylov", gcrNkrylov);
    read(paramtop, "VerboseP", verboseInner);
    read(paramtop, "InvTypeSloppy", invTypeSloppy);
  };

  void read(XMLReader& xml, const std::string& path, 
	    GCRInnerSolverParams& p)
  {
    GCRInnerSolverParams tmp(xml, path);
    p = tmp;
  }

  
  
  void write(XMLWriter& xml, const std::string& path, 
	     const GCRInnerSolverParams& p) {
    push(xml, path);
    write(xml, "RsdSloppy", p.tolSloppy);
    write(xml, "MaxIterSloppy", p.maxIterSloppy);
    write(xml, "NKrylov", p.gcrNkrylov);
    write(xml, "VerboseP", p.verboseInner);
    write(xml, "InvTypeSloppy", p.invTypeSloppy);
    pop(xml);

  }
}
