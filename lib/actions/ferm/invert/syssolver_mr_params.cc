// $Id: syssolver_mr_params.cc,v 1.2 2007-05-01 14:39:13 bjoo Exp $
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

    read(paramtop, "RsdMR", param.RsdMR);
    read(paramtop, "MaxMR", param.MaxMR);

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
    write(xml, "RsdMR", param.RsdMR);
    write(xml, "MaxMR", param.MaxMR);
    write(xml, "MROver", param.MROver);
    pop(xml);
  }

  //! Default constructor
  SysSolverMRParams::SysSolverMRParams()
  {
    RsdMR = zero;
    MaxMR = 0;
    MROver = 1.0;
  }

  //! Read parameters
  SysSolverMRParams::SysSolverMRParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
