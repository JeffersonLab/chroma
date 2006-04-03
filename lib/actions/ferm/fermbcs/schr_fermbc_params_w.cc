// $Id: schr_fermbc_params_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Schroedinger functional ferm boundary conditions
 */

#include "actions/ferm/fermbcs/schr_fermbc_params_w.h"

namespace Chroma 
{

  SchrFermBCParams::SchrFermBCParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SchrPhiMult", SchrPhiMult);
    read(paramtop, "decay_dir", decay_dir);
    read(paramtop, "loop_extent", loop_extent);
    read(paramtop, "theta", theta);
  }

  void read(XMLReader& xml, const std::string& path, SchrFermBCParams& p)   
  {
    SchrFermBCParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const SchrFermBCParams& p) 
  {
    push(xml, path);
    write(xml, "SchrPhiMult", p.SchrPhiMult);
    write(xml, "decay_dir", p.decay_dir);
    write(xml, "loop_extent", p.loop_extent);
    write(xml, "theta", p.theta);
    pop(xml);
  }

}
