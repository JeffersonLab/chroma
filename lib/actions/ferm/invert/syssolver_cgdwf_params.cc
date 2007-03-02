// $Id: syssolver_cgdwf_params.cc,v 1.1 2007-03-02 20:59:34 bjoo Exp $
/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_cgdwf_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverCGDWFParams& param)
  {
    XMLReader paramtop(xml, path);

    if( paramtop.count("RsdCGSingle") == 1 ) { 
      // Backward compatibility mode
      read(paramtop, "RsdCGSingle", param.RsdCGSingle);
    }
    else {
      read(paramtop, "RsdCG", param.RsdCGSingle);
    }

    if( paramtop.count("MaxCGSingle") == 1 ) { 
      // Backward compatibility mode
      read(paramtop, "MaxCGSingle", param.MaxCGSingle);
    }
    else {
      read(paramtop, "MaxCG", param.MaxCGSingle);
    }

    if( paramtop.count("RsdCGDouble") == 1 ) { 
      read(paramtop, "RsdCGDouble", param.RsdCGDouble);
    }
    else { 
      param.RsdCGDouble=param.RsdCGSingle;
    }

    if( paramtop.count("MaxCGDouble") == 1 ) { 
      read(paramtop, "MaxCGDouble", param.MaxCGDouble);
    }
    else { 
      param.MaxCGDouble=param.MaxCGSingle;
    }
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverCGDWFParams& param)
  {
    push(xml, path);

    write(xml, "invType", "CGDWF_INVERTER");
    write(xml, "RsdCGSingle", param.RsdCGSingle);
    write(xml, "MaxCGSingle", param.MaxCGSingle);
    write(xml, "RsdCGDouble", param.RsdCGDouble);
    write(xml, "MaxCGDouble", param.MaxCGDouble);


    pop(xml);
  }

  //! Default constructor
  SysSolverCGDWFParams::SysSolverCGDWFParams()
  {
    RsdCGSingle = RsdCGDouble = zero;
    MaxCGSingle = MaxCGDouble = 0;
  }

  //! Read parameters
  SysSolverCGDWFParams::SysSolverCGDWFParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
