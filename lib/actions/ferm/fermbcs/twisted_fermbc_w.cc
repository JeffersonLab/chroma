// $Id: twisted_fermbc_w.cc,v 1.1 2005-09-06 10:59:36 bjoo Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/twisted_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  // Name for both 4d and 5d 
  namespace TwistedFermBCEnv {
    const std::string name = "TWISTED_FERMBC";
  }

  // Readers and writerss for the params 
 //! Read parameters
  TwistedFermBCParams::TwistedFermBCParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "phases", boundary_phases);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, TwistedFermBCParams& param)
  {
    TwistedFermBCParams tmp(xml, path);
    param = tmp;
  }

  //! Write parameters
  void write(XMLWriter& xml_out, const string& path, const TwistedFermBCParams& param)
  {
    if ( path != "." )
      push(xml_out, path);
  
    write(xml_out, "FermBC", TwistedFermBCEnv::name);
    write(xml_out, "phases", param.boundary_phases);

    if( path != "." )
      pop(xml_out);
  }

  //! Name and registration
  namespace WilsonTypeTwistedFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      TwistedFermBCParams bc(xml_in, path);
      return new TwistedFermBC<LatticeFermion>(bc.boundary_phases);
    }

    //! Name to be used
    const std::string name = TwistedFermBCEnv::name;

    //! Register the fermbc
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
  }


  //! Name and registration
  namespace WilsonTypeTwistedFermBCArrayEnv
  {
    //! Callback function
    FermBC< multi1d<LatticeFermion> >* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      TwistedFermBCParams bc(xml_in, path);
      return new TwistedFermBC< multi1d<LatticeFermion> >(bc.boundary_phases);
    }

    //! Name to be used
    const std::string name = TwistedFermBCEnv::name;

    //! Register the fermbc
    const bool registered = TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBC);
  }

}
