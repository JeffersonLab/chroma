#ifndef __SYSSOLVER_QUDA_WILSON_PARAMS_H__
#define __SYSSOLVER_QUDA_WILSON_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include <string>
using namespace std;

namespace Chroma 
{
  struct SysSolverQUDAWilsonParams { 
    SysSolverQUDAWilsonParams(XMLReader& xml, const std::string& path);
    SysSolverQUDAWilsonParams() {
      solverType="CG";
    };
    SysSolverQUDAWilsonParams( const SysSolverQUDAWilsonParams& p) {
      WilsonParams = p.WilsonParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
      solverType = p.solverType;
      verboseP = p.verboseP;
    }
    WilsonFermActParams WilsonParams;
    bool AntiPeriodicT;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
    std::string solverType;
    bool verboseP;
  };

  void read(XMLReader& xml, const std::string& path, SysSolverQUDAWilsonParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverQUDAWilsonParams& param);



}

#endif


