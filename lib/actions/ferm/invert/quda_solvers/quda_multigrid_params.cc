#include "chromabase.h"
#include <string>
#include "actions/ferm/invert/quda_solvers/quda_multigrid_params.h"

using namespace QDP;

namespace Chroma {
  MULTIGRIDSolverParams::MULTIGRIDSolverParams(XMLReader& xml, 
					     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "Residual", tol);
    read(paramtop, "MaxIterations", maxIterations);
    read(paramtop, "SmootherType", smootherType);
    read(paramtop, "Verbosity", verbosity);
    read(paramtop, "Precision", prec);
    read(paramtop, "Reconstruct", reconstruct);
    read(paramtop, "NullVectors", nvec);
    read(paramtop, "Blocking", blocking);
    mg_levels = blocking.size()+1;
    read(paramtop, "GenerateNullspace", generate_nullspace);
    read(paramtop, "GenerateAllLevels", generate_all_levels);
    read(paramtop, "Pre-SmootherApplications", nu_pre);
    read(paramtop, "Post-SmootherApplications", nu_post);

    //FIXME: Should be an enum
    read(paramtop, "CycleType", cycle_type);
  

    // Read optional params otherwise use defaults
    if( paramtop.count("SchwarzType") > 0 ) { 
      read(paramtop, "SchwarzType", schwarzType);
    }

    // Read optional relaxation param
    if( paramtop.count("RelaxationOmegaMG") > 0 ) {
      read(paramtop, "RelaxationOmegaMG", relaxationOmegaMG);
    }
    else { 
	relaxationOmegaMG = Real(1.0);
    }

   // Read optional relaxation param
    if( paramtop.count("RelaxationOmegaOuter") > 0 ) {
      read(paramtop, "RelaxationOmegaOuter", relaxationOmegaOuter);
    }
    else {
        relaxationOmegaOuter = Real(1.0);
    }

    
  };

  void read(XMLReader& xml, const std::string& path, 
	    MULTIGRIDSolverParams& p)
  {
    MULTIGRIDSolverParams tmp(xml, path);
    p = tmp;
  }

  
  
  void write(XMLWriter& xml, const std::string& path, 
	     const MULTIGRIDSolverParams& p) {
    push(xml, path);
    write(xml, "Residual", p.tol);
    write(xml, "MaxIterations", p.maxIterations);
    write(xml, "SmootherType", p.smootherType);
    write(xml, "RelaxationOmegaMG", p.relaxationOmegaMG);
    write(xml, "RelaxationOmegaOuter", p.relaxationOmegaOuter);
    write(xml, "Verbosity", p.verbosity);
    write(xml, "Precision", p.prec);
    write(xml, "Reconstruct", p.reconstruct);
    write(xml, "SchwarzType", p.schwarzType);
    write(xml, "NullVectors", p.nvec);
    write(xml, "MultiGridLevels", p.mg_levels);
    write(xml, "GenerateNullSpace", p.generate_nullspace);
    write(xml, "GenerateAllLevels", p.generate_all_levels);
    write(xml, "CycleType", p.cycle_type);
    write(xml, "Pre-SmootherApplications", p.nu_pre);
    write(xml, "Post-SmootherApplications", p.nu_post);
    write(xml, "Blocking", p.blocking);
    pop(xml);

  }
}
