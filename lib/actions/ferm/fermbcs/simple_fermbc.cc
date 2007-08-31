 // $Id: simple_fermbc.cc,v 3.1 2007-08-31 03:33:46 edwards Exp $
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

    if (boundary.size() != Nd)
    {
      QDPIO::cerr << __func__ << ": boundary not the expected length = Nd" << endl;
      QDP_abort(1);
    }
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
