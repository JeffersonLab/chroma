#ifndef __MULTI_SYSSOLVER_QPHIX_CLOVER_PARAMS_H__
#define __MULTI_SYSSOLVER_QPHIX_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "io/enum_io/enum_type_map.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include <string>
#include "handle.h"

namespace Chroma 
{
  struct MultiSysSolverQPhiXCloverParams { 
    MultiSysSolverQPhiXCloverParams(XMLReader& xml, const std::string& path);
    MultiSysSolverQPhiXCloverParams() {
       RsdToleranceFactor = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
      VerboseP = false;
      SolutionCheckP = true;
      Delta = Real(-1);
    };

    MultiSysSolverQPhiXCloverParams( const MultiSysSolverQPhiXCloverParams& p) {
      CloverParams = p.CloverParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      Delta = p.Delta;
      MaxShifts = p.MaxShifts;
      RsdTarget = p.RsdTarget;
      VerboseP = p.VerboseP;
      SolutionCheckP = p.SolutionCheckP;
      RsdToleranceFactor = p.RsdToleranceFactor;
    }

   
    CloverFermActParams CloverParams;
    bool AntiPeriodicT;
    int MaxIter;
    int MaxShifts;
    multi1d<Real> RsdTarget;
    Real Delta;
    bool VerboseP;
    bool SolutionCheckP;
    Real RsdToleranceFactor;
  };

  void read(XMLReader& xml, const std::string& path, MultiSysSolverQPhiXCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const MultiSysSolverQPhiXCloverParams& param);



}

#endif


