/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/multi_syssolver_cg_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, MultiSysSolverCGParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);
  }

  // Writer parameters
  void write(XMLWriter& xml, const std::string& path, const MultiSysSolverCGParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "CG_INVERTER");
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);

    pop(xml);
  }

  //! Default constructor
  MultiSysSolverCGParams::MultiSysSolverCGParams()
  {
    RsdCG = zero;
    MaxCG = 0;
  }

  //! Read parameters
  MultiSysSolverCGParams::MultiSysSolverCGParams(XMLReader& xml, const std::string& path)
  {
    read(xml, path, *this);
  }

}
