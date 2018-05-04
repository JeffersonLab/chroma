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

    read(paramtop, "Reconstruct", reconstruct);
    
    read(paramtop, "Blocking", blocking);
    mg_levels = blocking.size()+1;
    nvec.resize(mg_levels-1);
    nu_pre.resize(mg_levels-1);
    nu_post.resize(mg_levels-1);

    read(paramtop, "NullVectors", nvec);
    read(paramtop, "Pre-SmootherApplications", nu_pre);
    read(paramtop, "Post-SmootherApplications", nu_post);
    if (nvec.size() != mg_levels-1 ) {
 
      QDPIO::cout<<"Warning. There are "<< blocking.size() 
		 << " blockings but only " << nvec.size() << " sets of NullVectors" << std::endl;
      QDP_abort(1);
    }

    if (nu_pre.size() != mg_levels-1 ) {
 
      QDPIO::cout<<"Error. There are "<< (mg_levels-1)  
		 << " blockings but only " << nu_pre.size() << " sets pre-smoothing iterations" << std::endl;
      QDP_abort(1);
    }

    if( paramtop.count("./OuterGCRNKrylov") == 1 ) {
	read(paramtop, "OuterGCRNKrylov", outer_gcr_nkrylov);
    }
    else { 
	outer_gcr_nkrylov = 12;
    }

    if( paramtop.count("./PrecondGCRNKrylov") == 1 ) { 
	read(paramtop, "PrecondGCRNKrylov", precond_gcr_nkrylov);
    }
    else {
        precond_gcr_nkrylov = 12;
    }

    if (nu_post.size() != mg_levels-1 ) {
 
      QDPIO::cout<<"Warning. There are "<< blocking.size() 
		 << " blockings but only " << nu_post.size() << " sets post-smoothing iterations " << std::endl;
      QDP_abort(1);
    }

    read(paramtop, "GenerateNullspace", generate_nullspace);
    read(paramtop, "GenerateAllLevels", generate_all_levels);

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
    /* FIXME: This should go in the general solver interface, and work for all GCR solvers, not just GCR inner params */
    write(xml, "OuterGCRNKrylov", p.outer_gcr_nkrylov);
    write(xml, "PrecondGCRNKrylov", p.precond_gcr_nkrylov);
    write(xml, "Blocking", p.blocking);
    pop(xml);

  }
}
