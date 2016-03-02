#ifndef __QUDA_MULTIGRID_PARAMS_H__
#define __QUDA_MULTIGRID_PARAMS_H__

#include "chromabase.h"
#include <string>
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{

  struct MULTIGRIDSolverParams {
    
    Real tol;
    int  maxIterations;
    QudaSolverType smootherType;
    bool verbosity;
    QudaPrecisionType prec;
    QudaReconsType reconstruct;
    QudaSchwarzMethod schwarzType;
    multi1d<int> nvec;
    int mg_levels;
    bool generate_nullspace;
    bool generate_all_levels;
    multi1d<int> nu_pre;
    multi1d<int> nu_post;
    multi1d< multi1d<int> > blocking;
    std::string cycle_type;
    Real relaxationOmegaMG;
    Real relaxationOmegaOuter;
 
    MULTIGRIDSolverParams(XMLReader& xml, const std::string& path);
    MULTIGRIDSolverParams() {
      tol = .000001;
      relaxationOmegaMG =Real(1.0);
      relaxationOmegaOuter = Real(1.0);
      maxIterations = 10;
      smootherType = MR;
      verbosity = false;
      prec = DEFAULT;
      reconstruct = RECONS_NONE;
      schwarzType = ADDITIVE_SCHWARZ;

      mg_levels = 2;
      blocking.resize(mg_levels-1);
      nvec.resize(mg_levels-1);
      nu_pre.resize(mg_levels-1);
      nu_post.resize(mg_levels-1);
      generate_nullspace = true;
      for(int l = 0; l < mg_levels - 1; l++) 
      {
	blocking[l].resize(4);
        blocking[l][0] = blocking[l][1] = blocking[l][2] = blocking[l][3] = 4;
	nu_pre[l] = 2;
	nu_post[l] = 2;
	nvec[l] = 16;
      }
      generate_all_levels = true;
      cycle_type = "MG_VCYCLE";

    };

  };
  void read(XMLReader& xml, const std::string& path, MULTIGRIDSolverParams& p);
 
  void write(XMLWriter& xml, const std::string& path, 
	     const MULTIGRIDSolverParams& param);

}
#endif
