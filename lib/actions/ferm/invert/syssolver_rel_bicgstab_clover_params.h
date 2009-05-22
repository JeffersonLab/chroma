#ifndef __SYSSOLVER_REL_BICGSTAB_PARAMS_H__
#define __SYSSOLVER_REL_BICGSTAB_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{
  struct SysSolverReliableBiCGStabCloverParams { 
    SysSolverReliableBiCGStabCloverParams(XMLReader& xml, const std::string& path);
    SysSolverReliableBiCGStabCloverParams() {};
    SysSolverReliableBiCGStabCloverParams( const SysSolverReliableBiCGStabCloverParams& p) {
      clovParams = p.clovParams;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
    }
    CloverFermActParams clovParams;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
  };

  typedef SysSolverReliableBiCGStabCloverParams SysSolverReliableCGCloverParams;
  void read(XMLReader& xml, const std::string& path, SysSolverReliableBiCGStabCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverReliableBiCGStabCloverParams& param);



}

#endif


