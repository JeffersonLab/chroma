// $Id: simple_fermbc.cc,v 1.1 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Name
  namespace SimpleFermBCEnv
  {
    const std::string name = "SIMPLE_FERMBC";
  };

  //! Read parameters
  SimpleFermBCParams::SimpleFermBCParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "boundary", boundary);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, SimpleFermBCParams& param)
  {
    SimpleFermBCParams tmp(xml, path);
    param = tmp;
  }

  //! Write parameters
  void write(XMLWriter& xml_out, const string& path, const SimpleFermBCParams& param)
  {
    if ( path != "." )
      push(xml_out, path);
  
    write(xml_out, "FermBC", SimpleFermBCEnv::name);
    write(xml_out, "boundary", param.boundary);

    if( path != "." )
      pop(xml_out);
  }

}
