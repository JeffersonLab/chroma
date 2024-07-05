#ifndef __QUDA_MULTIGRID_PARAMS_H__
#define __QUDA_MULTIGRID_PARAMS_H__

#include "chromabase.h"
#include <string>
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{

  struct MULTIGRIDSolverParams {
    
	// tol, max iterations, and smootherType are per level.
    multi1d<Real> tol;
    multi1d<int> maxIterations;
    multi1d<QudaSolverType> smootherType;
    multi1d<QudaSolverType> coarseSolverType;
    multi1d<Real> smootherTol;
    multi1d<Real> relaxationOmegaMG;
    multi1d<QudaSchwarzMethod> smootherSchwarzType;
    multi1d<int> smootherSchwarzCycle;
    multi1d<QudaPrecisionType> smootherHaloPrecision;

    bool verbosity;
    QudaPrecisionType prec;
    QudaReconsType reconstruct;
    QudaSchwarzMethod schwarzType;
    multi1d<int> nvec;
    int mg_levels; 
    bool generate_nullspace;
    bool generate_all_levels;
    bool check_multigrid_setup;
    multi1d<bool> setup_on_gpu;
    multi1d<int> nu_pre;
    multi1d<int> nu_post;
    multi1d< multi1d<int> > blocking;
		multi1d<int> nvec_batch;
    int outer_gcr_nkrylov;
    int precond_gcr_nkrylov;
    std::string cycle_type;
    Real relaxationOmegaOuter;
    multi1d<QudaSolverType> subspaceSolver;
    multi1d<int> maxIterSubspaceCreate;
    multi1d<Real> rsdTargetSubspaceCreate;
    multi1d<int> maxIterSubspaceRefresh;
    
    MULTIGRIDSolverParams(XMLReader& xml, const std::string& path);
    MULTIGRIDSolverParams() {

      relaxationOmegaMG =Real(1.0);
      relaxationOmegaOuter = Real(1.0);
      maxIterations = 10;
      smootherType = MR;
      smootherHaloPrecision = DEFAULT;
      verbosity = false;
      prec = DEFAULT;
      reconstruct = RECONS_NONE;
      schwarzType = ADDITIVE_SCHWARZ;

      mg_levels = 2;
      setup_on_gpu.resize(mg_levels-1);
      // Default is to set up on gpu
      for(int i=0; i < mg_levels - 1; ++i) { 
        setup_on_gpu[i] = true;
      }
      blocking.resize(mg_levels-1);
      nvec.resize(mg_levels-1);
			nvec_batch.resize(mg_levels-1);
      nu_pre.resize(mg_levels-1);
      nu_post.resize(mg_levels-1);
      maxIterSubspaceCreate.resize(mg_levels-1);
      maxIterSubspaceRefresh.resize(mg_levels-1);
      rsdTargetSubspaceCreate.resize(mg_levels-1);
      tol.resize(mg_levels);
      maxIterations.resize(mg_levels);
      coarseSolverType.resize(mg_levels);

      smootherType.resize(mg_levels);
      smootherTol.resize(mg_levels);
      relaxationOmegaMG.resize(mg_levels);
      smootherSchwarzType.resize(mg_levels);
      smootherSchwarzCycle.resize(mg_levels);
      generate_nullspace = true;
      for(int l = 0; l < mg_levels - 1; l++) 
      {
    	  blocking[l].resize(4);
    	  blocking[l][0] = blocking[l][1] = blocking[l][2] = blocking[l][3] = 4;
    	  nu_pre[l] = 2;
    	  nu_post[l] = 2;
    	  nvec[l] = 16;
				nvec_batch[l]=1; // the batch size for Nvec solves is 1 by default
    	  // Default params:
    	  maxIterSubspaceCreate[l] = 500;
    	  rsdTargetSubspaceCreate[l] = 5.0e-6;
    	  maxIterSubspaceRefresh[l] = maxIterSubspaceCreate[l];

    	  tol[l] = 1.0e-5;
    	  coarseSolverType[l] = GCR;
    	  maxIterations[l] = 12;
      }
      for(int l=0; l < mg_levels; ++l) {
    	  smootherType[l] = MR;
    	  smootherTol[l] = Real(0.25);
    	  relaxationOmegaMG[l] = 0.85;
    	  smootherSchwarzType[l] = INVALID_SCHWARZ;
    	  smootherSchwarzCycle[l] = 1;
      }
      outer_gcr_nkrylov = 12;
      precond_gcr_nkrylov = 12;

      generate_all_levels = true;
      check_multigrid_setup= true;
      cycle_type = "MG_VCYCLE";

    };

  };
  void read(XMLReader& xml, const std::string& path, MULTIGRIDSolverParams& p);
 
  void write(XMLWriter& xml, const std::string& path, 
	     const MULTIGRIDSolverParams& param);

}
#endif
