// $Id: schr_gaugebc_params.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Schroedinger functional gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{

  SchrGaugeBCParams::SchrGaugeBCParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SchrFunType", SchrFun);
    read(paramtop, "SchrPhiMult", SchrPhiMult);
  }

  void read(XMLReader& xml, const std::string& path, SchrGaugeBCParams& p)   
  {
    SchrGaugeBCParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const SchrGaugeBCParams& p) 
  {
    push(xml, path);
    write(xml, "SchrFunType", p.SchrFun);
    write(xml, "SchrPhiMult", p.SchrPhiMult);
    pop(xml);
  }

}
