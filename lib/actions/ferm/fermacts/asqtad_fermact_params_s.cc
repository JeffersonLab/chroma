// $Id: asqtad_fermact_params_s.cc,v 3.1 2006-11-16 19:49:33 kostas Exp $
/*! \file
 *  \brief Asqtad fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/asqtad_fermact_params_s.h"

#include "io/param_io.h"

namespace Chroma
{
  //! Default constructor
  AsqtadFermActParams::AsqtadFermActParams()
  {
    Mass = 0.0;
    u0  = 1.0;
  }


  //! Read parameters
  AsqtadFermActParams::AsqtadFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "Mass", Mass);
    read(paramtop, "u0", u0);
    
    //  Read optional anisoParam.
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, AsqtadFermActParams& param)
  {
    AsqtadFermActParams tmp(xml, path);
    param = tmp;
  }

  //! Writer parameters
  void write(XMLWriter& xml, const string& path, const AsqtadFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    write(xml, "u0", param.u0);
    pop(xml);
  }
}
