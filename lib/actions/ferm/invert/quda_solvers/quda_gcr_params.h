#ifndef __QUDA_GCR_PARAMS_H__
#define __QUDA_GCR_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{

  struct GCRInnerSolverParams {
    
    Real tolPrecondition;
    int  maxIterPrecondition;
    int  gcrNkrylov;
    bool verboseInner;
    QudaSolverType invTypePrecondition;
    QudaPrecisionType precPrecondition;
    QudaReconsType reconstructPrecondition;
    QudaSchwarzMethod schwarzType;
    int preconditionCycle;
    
    GCRInnerSolverParams(XMLReader& xml, const std::string& path);
    GCRInnerSolverParams() {
      tolPrecondition=0;
      maxIterPrecondition=0;
      gcrNkrylov=0;
      verboseInner=false;
      invTypePrecondition=MR;
      precPrecondition=DEFAULT;
      reconstructPrecondition=RECONS_NONE;
      schwarzType = ADDITIVE_SCHWARZ;
      preconditionCycle = 1;
    };

  };
  void read(XMLReader& xml, const std::string& path, GCRInnerSolverParams& p);
 
  void write(XMLWriter& xml, const std::string& path, 
	     const GCRInnerSolverParams& param);

}
#endif
