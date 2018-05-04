/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_mr_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, SysSolverMRParams& param)
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
  void write(XMLWriter& xml, const std::string& path, const SysSolverMRParams& param)
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
  SysSolverMRParams::SysSolverMRParams(XMLReader& xml, const std::string& path)
  {
    read(xml, path, *this);
  }

}
