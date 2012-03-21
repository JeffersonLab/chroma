#ifndef __QUDA_GCR_PARAMS_H__
#define __QUDA_GCR_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{

  struct GCRInnerSolverParams {

    Real tolSloppy;
    int  maxIterSloppy;
    int  gcrNkrylov;
    bool verboseInner;
    QudaSolverType invTypeSloppy;

    GCRInnerSolverParams(XMLReader& xml, const std::string& path);
    GCRInnerSolverParams() {
      tolSloppy=0;
      maxIterSloppy=0;
      gcrNkrylov=0;
      verboseInner=false;
      invTypeSloppy=CG;
    };

  };
  void read(XMLReader& xml, const std::string& path, GCRInnerSolverParams& p);
 
  void write(XMLWriter& xml, const std::string& path, 
	     const GCRInnerSolverParams& param);

}
#endif
