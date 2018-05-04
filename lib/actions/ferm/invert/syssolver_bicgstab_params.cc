/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_bicgstab_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, SysSolverBiCGStabParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdBiCGStab", param.RsdBiCGStab);
    read(paramtop, "MaxBiCGStab", param.MaxBiCGStab);

  }

  // Writer parameters
  void write(XMLWriter& xml, const std::string& path, const SysSolverBiCGStabParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "BICGSTAB_INVERTER");
    write(xml, "RsdBiCGStab", param.RsdBiCGStab);
    write(xml, "MaxBiCGStab", param.MaxBiCGStab);
    pop(xml);
  }

  //! Default constructor
  SysSolverBiCGStabParams::SysSolverBiCGStabParams()
  {
    RsdBiCGStab = zero;
    MaxBiCGStab = 0;
  }

  //! Read parameters
  SysSolverBiCGStabParams::SysSolverBiCGStabParams(XMLReader& xml, const std::string& path)
  {
    read(xml, path, *this);
  }

}
