#ifndef __SYSSOLVER_CG_QUDA_WILSON_PARAMS_H__
#define __SYSSOLVER_CG_QUDA_WILSON_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
namespace Chroma 
{
  struct SysSolverCGQUDAWilsonParams { 
    SysSolverCGQUDAWilsonParams(XMLReader& xml, const std::string& path);
    SysSolverCGQUDAWilsonParams() {};
    SysSolverCGQUDAWilsonParams( const SysSolverCGQUDAWilsonParams& p) {
      WilsonParams = p.WilsonParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
    }
    WilsonFermActParams WilsonParams;
    bool AntiPeriodicT;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
  };

  void read(XMLReader& xml, const std::string& path, SysSolverCGQUDAWilsonParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverCGQUDAWilsonParams& param);



}

#endif


