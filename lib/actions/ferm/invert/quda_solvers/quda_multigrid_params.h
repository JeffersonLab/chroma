#ifndef __QUDA_MULTIGRID_PARAMS_H__
#define __QUDA_MULTIGRID_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{

  struct MULTIGRIDSolverParams {
    
    Real tol;
    int  maxIterations;
    bool verbosity;
    QudaSolverType invType;
    QudaPrecisionType prec;
    QudaReconsType reconstruct;
    QudaSchwarzMethod schwarzType;
    int nvec;
    int mg_levels;
    bool generate_nullspace;
    int nu_pre;
    int nu_post;
    multi1d< multi1d<int> > blocking;
    
    MULTIGRIDSolverParams(XMLReader& xml, const std::string& path);
    MULTIGRIDSolverParams() {
      tol = .000001;
      maxIterations = 10;
      verbosity = false;
      invType = CG;
      prec = DEFAULT;
      reconstruct = RECONS_NONE;
      schwarzType = ADDITIVE_SCHWARZ;
      nvec = 16;      
      mg_levels = 2;
      generate_nullspace = true;
      nu_pre = 2;
      nu_post = 2;
      for(int l = 0; l < mg_levels - 1; l++) 
      {
	blocking[l].resize(4);
        blocking[l][0] = blocking[l][1] = blocking[l][2] = blocking[l][3] = 4;
      }
    };

  };
  void read(XMLReader& xml, const std::string& path, MULTIGRIDSolverParams& p);
 
  void write(XMLWriter& xml, const std::string& path, 
	     const MULTIGRIDSolverParams& param);

}
#endif
