#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/quda_multigrid_params.h"

using namespace QDP;

namespace Chroma {
  MULTIGRIDSolverParams::MULTIGRIDSolverParams(XMLReader& xml, 
					     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "Residual", tol);
    read(paramtop, "MaxIterations", maxIterations);
    read(paramtop, "InverterType", invType);
    read(paramtop, "Precision", prec);
    read(paramtop, "Reconstruct", reconstruct);
    read(paramtop, "NullVectors", nvec);
    read(paramtop, "MultiGridLevels", mg_levels);
    read(paramtop, "GenerateNullspace", generate_nullspace);
    read(paramtop, "Pre-SmootherApplications", nu_pre);
    read(paramtop, "Post-SmootherApplications", nu_post);
    read(paramtop, "Blocking", blocking);

    // Read optional params otherwise use defaults
    if( paramtop.count("SchwarzType") > 0 ) { 
      read(paramtop, "SchwarzType", schwarzType);
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
    write(xml, "InverterType", p.invType);
    write(xml, "Precision", p.prec);
    write(xml, "Reconstruct", p.reconstruct);
    write(xml, "SchwarzType", p.schwarzType);
    write(xml, "NullVectors", p.nvec);
    write(xml, "MultiGridLevels", p.mg_levels);
    write(xml, "GenerateNullSpace", p.generate_nullspace);
    write(xml, "Pre-SmootherApplications", p.nu_pre);
    write(xml, "Post-SmootherApplications", p.nu_post);
    write(xml, "Blocking", p.blocking);
    pop(xml);

  }
}
