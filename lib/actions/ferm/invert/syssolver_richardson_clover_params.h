#ifndef __SYSSOLVER_RICHARDSON_PARAMS_H__
#define __SYSSOLVER_RICHARDSON_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{
  struct SysSolverRichardsonCloverParams { 
    SysSolverRichardsonCloverParams(XMLReader& xml, const std::string& path);
    SysSolverRichardsonCloverParams() {};
    SysSolverRichardsonCloverParams( const SysSolverRichardsonCloverParams& p) {
      clovParams = p.clovParams;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      innerSolverParams = p.innerSolverParams;
    }
    CloverFermActParams clovParams;
    int MaxIter;
    Real RsdTarget;
    GroupXML_t innerSolverParams;
  };


  void read(XMLReader& xml, const std::string& path, SysSolverRichardsonCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverRichardsonCloverParams& param);



}

#endif


