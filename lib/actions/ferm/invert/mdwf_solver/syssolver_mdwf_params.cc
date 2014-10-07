#include "actions/ferm/invert/mdwf_solver/syssolver_mdwf_params.h"
#include "chromabase.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {

  SysSolverMDWFParams::SysSolverMDWFParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    if( paramtop.count("MaxIterRestart") > 0 ) { 
      read(paramtop, "MaxIterRestart", MaxIterRestart);
    }
    else {
      MaxIterRestart = MaxIter;
    }


    read(paramtop, "RsdTarget", RsdTarget);
    if( paramtop.count("RsdTargetRestart") > 0 ) { 
      read(paramtop, "RsdTargetRestart", RsdTargetRestart);
    }
    else {
      RsdTargetRestart = RsdTarget;
    }

    read(paramtop, "N5", N5);
    read(paramtop, "Mass", Mass);
    read(paramtop, "OverMass", OverMass);
    if(  (paramtop.count("c5") == 0)  && (paramtop.count("b5") == 0) ) {
      // No b5 and no c5 -- use shamir defaults
      b5=Real(1); 
      c5=Real(0);
    }
    else { 
      read(paramtop, "b5", b5);
      read(paramtop, "c5", c5);
    }

    if ( paramtop.count("AnisoParam") > 0 ) { 
      read(paramtop, "AnisoParam", anisoParam);
    }
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverMDWFParams& p)
  {
    SysSolverMDWFParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverMDWFParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "MaxIterRestart", p.MaxIterRestart);
    write(xml, "RsdTargetRestart", p.RsdTargetRestart);

    write(xml, "N5", p.N5);
    write(xml, "Mass", p.Mass);
    write(xml, "OverMass", p.OverMass);

    write(xml, "b5", p.b5);
    write(xml, "c5", p.c5);
    write(xml, "AnisoParam", p.anisoParam);

    pop(xml);

  }



}
