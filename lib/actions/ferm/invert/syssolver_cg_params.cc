// $Id: syssolver_cg_params.cc,v 3.2 2006-10-15 04:17:00 edwards Exp $
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

    if (paramtop.count("numRestarts") > 0)
      read(paramtop, "numRestarts", param.numRestarts);
    else
      param.numRestarts = 1;

    if (param.numRestarts <= 0)
    {
      QDPIO::cerr << __func__ << ": invalid SysSolverCGParams::numRestarts = " 
		  << param.numRestarts << endl;
      QDP_abort(1);
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
    write(xml, "numRestarts", param.numRestarts);

    pop(xml);
  }

  //! Default constructor
  SysSolverCGParams::SysSolverCGParams()
  {
    RsdCG = zero;
    MaxCG = 0;
    numRestarts = 1;
  }

  //! Read parameters
  SysSolverCGParams::SysSolverCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
