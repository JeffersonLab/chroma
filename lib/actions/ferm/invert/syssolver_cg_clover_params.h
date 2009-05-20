#ifndef __SYSSOLVER_CG_CLOVER_PARAMS_H__
#define __SYSSOLVER_CG_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{
  struct SysSolverCGCloverParams { 
    SysSolverCGCloverParams(XMLReader& xml, const std::string& path);
    SysSolverCGCloverParams() {};
    SysSolverCGCloverParams( const SysSolverCGCloverParams& p) {
      clovParams = p.clovParams;
      MaxCG = p.MaxCG;
      RsdCG = p.RsdCG;
    }
    CloverFermActParams clovParams;
    int MaxCG;
    Real RsdCG;
  };


  void read(XMLReader& xml, const std::string& path, SysSolverCGCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverCGCloverParams& param);



}

#endif


