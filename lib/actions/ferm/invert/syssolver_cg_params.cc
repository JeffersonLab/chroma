// $Id: syssolver_cg_params.cc,v 3.5 2008-03-25 10:43:44 mcneile Exp $
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


    int aa = paramtop.count("MinCG") ;
    if( aa  > 0 ) { 
      read(paramtop,"MinCG",param.MinCG);
    }
    else {
      param.MinCG = 0 ; 
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
    write(xml, "MinCG", param.MinCG);
    write(xml, "RsdCGRestart", param.RsdCGRestart);
    write(xml, "MaxCGRestart", param.MaxCGRestart);
    pop(xml);
  }

  //! Default constructor
  SysSolverCGParams::SysSolverCGParams()
  {
    RsdCG = zero;
    MaxCG = 0;
    MinCG = 0;
    RsdCGRestart = RsdCG;
    MaxCGRestart = MaxCG;
  }

  //! Read parameters
  SysSolverCGParams::SysSolverCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
