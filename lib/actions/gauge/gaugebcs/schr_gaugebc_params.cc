// $Id: schr_gaugebc_params.cc,v 2.2 2006-03-13 05:19:01 edwards Exp $
/*! \file
 *  \brief Schroedinger functional gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{

  SchrGaugeBCParams::SchrGaugeBCParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SchrPhiMult", SchrPhiMult);
    read(paramtop, "decay_dir", decay_dir);
    read(paramtop, "loop_extent", loop_extent);
  }

  void read(XMLReader& xml, const std::string& path, SchrGaugeBCParams& p)   
  {
    SchrGaugeBCParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const SchrGaugeBCParams& p) 
  {
    push(xml, path);
    write(xml, "SchrPhiMult", p.SchrPhiMult);
    write(xml, "decay_dir", p.decay_dir);
    write(xml, "loop_extent", p.loop_extent);
    pop(xml);
  }

}
