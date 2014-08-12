#ifndef __SYSSOLVER_MDWF_PARAMS_H__
#define __SYSSOLVER_MDWF_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include <string>
#include "io/aniso_io.h"

#include "handle.h"
using namespace std;

namespace Chroma 
{
  struct SysSolverMDWFParams { 
    SysSolverMDWFParams(XMLReader& xml, const std::string& path);
    SysSolverMDWFParams() {
      // Default is shamir;
      b5=Real(1);
      c5=Real(0);
      anisoParam.anisoP=false;
    };

    SysSolverMDWFParams( const SysSolverMDWFParams& p) {
      MaxIter = p.MaxIter;
      MaxIterRestart = p.MaxIterRestart;
      RsdTarget = p.RsdTarget;
      RsdTargetRestart = p.RsdTargetRestart;

      N5 = p.N5;
      Mass = p.Mass;
      OverMass = p.OverMass;
      c5=p.c5;
      b5=p.b5;
      anisoParam = p.anisoParam;
    }

    int MaxIter;
    int MaxIterRestart;

    Real RsdTarget;
    Real RsdTargetRestart;
    int N5; 
    Real Mass;
    Real OverMass;
    Real c5;
    Real b5;
    AnisoParam_t anisoParam;
  };

  void read(XMLReader& xml, const std::string& path, SysSolverMDWFParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverMDWFParams& param);



}

#endif


