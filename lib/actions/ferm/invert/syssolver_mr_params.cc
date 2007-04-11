// $Id: syssolver_mr_params.cc,v 1.1 2007-04-11 03:42:07 edwards Exp $
/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_mr_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverMRParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);

    if (paramtop.count("MROver") > 0)
      read(paramtop, "MROver", param.MROver);
    else
      param.MROver = 1.0;
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverMRParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "MR_INVERTER");
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "MROver", param.MROver);
    pop(xml);
  }

  //! Default constructor
  SysSolverMRParams::SysSolverMRParams()
  {
    RsdCG = zero;
    MaxCG = 0;
    MROver = 1.0;
  }

  //! Read parameters
  SysSolverMRParams::SysSolverMRParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
