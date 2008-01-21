// $Id: wilson_coarse_fine_fermact_params_w.cc,v 3.1 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Wilson coarse-fine 2+2 anisotropic lattice fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/wilson_coarse_fine_fermact_params_w.h"

namespace Chroma
{
  //! Default constructor
  WilsonCoarseFineFermActParams::WilsonCoarseFineFermActParams()
  {
    Mass = 0.0;
    gamma_f = 1.0;
    coarse_dirs.resize(Nd);
    coarse_dirs = true;
  }


  //! Read parameters
  WilsonCoarseFineFermActParams::WilsonCoarseFineFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "Mass", Mass);
    read(paramtop, "coarse_dirs", coarse_dirs);
    read(paramtop, "gamma_f", gamma_f);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, WilsonCoarseFineFermActParams& param)
  {
    WilsonCoarseFineFermActParams tmp(xml, path);
    param = tmp;
  }

  //! Writer parameters
  void write(XMLWriter& xml, const string& path, const WilsonCoarseFineFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    write(xml, "coarse_dirs", param.coarse_dirs);
    write(xml, "gamma_f", param.gamma_f);

    pop(xml);
  }
}
