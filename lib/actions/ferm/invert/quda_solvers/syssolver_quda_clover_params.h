#ifndef __SYSSOLVER_QUDA_CLOVER_PARAMS_H__
#define __SYSSOLVER_QUDA_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
namespace Chroma 
{
  struct SysSolverQUDACloverParams { 
    SysSolverQUDACloverParams(XMLReader& xml, const std::string& path);
    SysSolverQUDACloverParams() {};
    SysSolverQUDACloverParams( const SysSolverQUDACloverParams& p) {
      CloverParams = p.CloverParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
    }
    CloverFermActParams CloverParams;
    bool AntiPeriodicT;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
  };

  void read(XMLReader& xml, const std::string& path, SysSolverQUDACloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverQUDACloverParams& param);



}

#endif


