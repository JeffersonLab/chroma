#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/quda_gcr_params.h"

using namespace QDP;

namespace Chroma {
  GCRInnerSolverParams::GCRInnerSolverParams(XMLReader& xml, 
					     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "RsdPrecondition", tolPrecondition);
    read(paramtop, "MaxIterPrecondition", maxIterPrecondition);
    read(paramtop, "NKrylov", gcrNkrylov);
    read(paramtop, "VerboseP", verboseInner);
    read(paramtop, "InvTypePrecondition", invTypePrecondition);
    read(paramtop, "PrecPrecondition", precPrecondition);
    read(paramtop, "ReconstructPrecondition", reconstructPrecondition);

    // Read GCR Specific params
    // Read optional GCR params otherwise use defaults
    if( paramtop.count("SchwarzType") > 0 ) { 
      read(paramtop, "SchwarzType", schwarzType);
    }
    
    if( paramtop.count("PreconditionCycle") > 0 ) { 
      read(paramtop, "PreconditionCycle", preconditionCycle);
    }
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
    write(xml, "RsdPrecondition", p.tolPrecondition);
    write(xml, "MaxIterPrecondition", p.maxIterPrecondition);
    write(xml, "NKrylov", p.gcrNkrylov);
    write(xml, "VerboseP", p.verboseInner);
    write(xml, "InvTypePrecondition", p.invTypePrecondition);
    write(xml, "SchwarzType", p.schwarzType);
    write(xml, "PreconditionCycle", p.preconditionCycle);
    write(xml, "PrecPrecondition", p.precPrecondition);
    write(xml, "ReconstructPrecondition", p.reconstructPrecondition);
    pop(xml);

  }
}
