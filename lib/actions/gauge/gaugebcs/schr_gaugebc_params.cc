// $Id: schr_gaugebc_params.cc,v 3.1 2006-04-11 17:24:59 edwards Exp $
/*! \file
 *  \brief Schroedinger functional gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{

  SchrGaugeBCParams::SchrGaugeBCParams()
  {
    SchrPhiMult = 1;
    decay_dir = Nd-1;
    loop_extent = 1;
  }

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
