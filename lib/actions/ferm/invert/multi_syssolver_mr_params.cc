/*! \file
 *  \brief Params of MR inverter
 */

#include "actions/ferm/invert/multi_syssolver_mr_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, MultiSysSolverMRParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);
  }

  // Writer parameters
  void write(XMLWriter& xml, const std::string& path, const MultiSysSolverMRParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "MR_INVERTER");
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);

    pop(xml);
  }

  //! Default constructor
  MultiSysSolverMRParams::MultiSysSolverMRParams()
  {
    RsdCG = zero;
    MaxCG = 0;
  }

  //! Read parameters
  MultiSysSolverMRParams::MultiSysSolverMRParams(XMLReader& xml, const std::string& path)
  {
    read(xml, path, *this);
  }

}
