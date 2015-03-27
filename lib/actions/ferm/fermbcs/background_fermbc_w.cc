// $Id: twisted_fermbc_w.cc,v 3.1 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/background_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  // Readers and writerss for the params 
  //! Read parameters
  BackgroundFermBCParams::BackgroundFermBCParams(XMLReader& xml, const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // The BASE Simple BC boundary
    read(paramtop, "boundary", boundary);
    read(paramtop, "gamma", gamma);
    read(paramtop, "lambda", lambda);

    if( boundary.size() != Nd ) { 
      QDPIO::cerr << "BackgroundFermBCParams: Invalid size for boundary. Should be " << Nd << "  but is " << boundary.size() << std::endl;
      QDP_abort(1);
    }

    if(( gamma < 0 )||(gamma>15)) { 
      QDPIO::cerr << "BackgroundFermBCParams: gamma should be [0,15]" ;
      QDPIO::cerr<< "  but is " << gamma << std::endl;
      QDP_abort(1);
    }
  }

  //! Read parameters
  void read(XMLReader& xml, const std::string& path, BackgroundFermBCParams& param)
  {
    BackgroundFermBCParams tmp(xml, path);
    param = tmp;
  }

  //! Write parameters
  void write(XMLWriter& xml_out, const std::string& path, const BackgroundFermBCParams& param)
  {
    if ( path != "." )
      push(xml_out, path);
  
    write(xml_out, "FermBC", WilsonTypeBackgroundFermBCEnv::name);
    write(xml_out, "boundary", param.boundary);
    write(xml_out, "gamma", param.gamma);
    write(xml_out, "lambda", param.lambda);

    if( path != "." )
      pop(xml_out);
  }


  //! Name and registration
  namespace WilsonTypeBackgroundFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      BackgroundFermBCParams bc(xml_in, path);
      return new BackgroundFermBC<LatticeFermion>(bc.gamma,
						  bc.lambda,
						  bc.boundary);
    }

    //! Name to be used
    const std::string name = "BACKGROUND_FERMBC";

    static bool registered = false;

    //! Register all the factories
    // Register all objects
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
	registered = true;
      }
      return success;
    }
  }
}
