// $Id: syssolver_eigcg_params.cc,v 1.1 2007-09-25 21:17:12 edwards Exp $
/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverCGParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);
    
    if( paramtop.count("RsdCGRestart") > 0 ) { 
      read(paramtop, "RsdCGRestart", param.RsdCGRestart);
    }
    else {
      param.RsdCGRestart = param.RsdCG;
    }

    if( paramtop.count("MaxCGRestart") > 0 ) { 
      read(paramtop, "MaxCGRestart", param.MaxCGRestart);
    }
    else {
      param.MaxCGRestart = param.MaxCG;
    }
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverCGParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "CG_INVERTER");
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "RsdCGRestart", param.RsdCGRestart);
    write(xml, "MaxCGRestart", param.RsdCGRestart);
    pop(xml);
  }

  //! Default constructor
  SysSolverCGParams::SysSolverCGParams()
  {
    RsdCG = zero;
    MaxCG = 0;
    RsdCGRestart = RsdCG;
    MaxCGRestart = MaxCG;
  }

  //! Read parameters
  SysSolverCGParams::SysSolverCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
