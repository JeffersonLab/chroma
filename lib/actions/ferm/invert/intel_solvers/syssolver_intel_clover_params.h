#ifndef __SYSSOLVER_INTEL_CLOVER_PARAMS_H__
#define __SYSSOLVER_INTEL_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include <string>
#include "handle.h"
using namespace std;

namespace Chroma 
{
  struct SysSolverIntelCloverParams { 
    SysSolverIntelCloverParams(XMLReader& xml, const std::string& path);
    SysSolverIntelCloverParams() {
       RsdToleranceFactor = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
      TuneP = false ; //< v0.3 autotune feature
      VerboseP = false;
      MinCt = 1;
      CompressP = false;
    };

    SysSolverIntelCloverParams( const SysSolverIntelCloverParams& p) {
      CloverParams = p.CloverParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      VerboseP = p.VerboseP;
      NCores = p.NCores;
      ThreadsPerCore = p.ThreadsPerCore;
      By=p.By;
      Bz=p.Bz;
      Sy=p.Sy;
      Sz=p.Sz;
      PadXY = p.PadXY;
      PadXYZ = p.PadXYZ;
      RsdToleranceFactor = p.RsdToleranceFactor;
      TuneP = p.TuneP;
      CompressP = p.CompressP;
      MinCt = p.MinCt;
    }

   
    CloverFermActParams CloverParams;
    bool AntiPeriodicT;
    int MaxIter;
    Real RsdTarget;
    bool VerboseP;
    int NCores;
    int ThreadsPerCore;
    int By;
    int Bz;
    int Sy;
    int Sz;
    int PadXY;
    int PadXYZ;
    int MinCt;
    Real RsdToleranceFactor;
    bool TuneP;
    bool CompressP;
  };

  void read(XMLReader& xml, const std::string& path, SysSolverIntelCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverIntelCloverParams& param);



}

#endif


