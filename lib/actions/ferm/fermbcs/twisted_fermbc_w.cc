// $Id: twisted_fermbc_w.cc,v 2.0 2005-09-25 21:04:27 edwards Exp $
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

    // The BASE Simple BC boundary
    read(paramtop, "boundary", boundary);
    read(paramtop, "phases_by_pi", phases_by_pi);
    read(paramtop, "phases_dir", phases_dir);

    if( boundary.size() != Nd ) { 
      QDPIO::cerr << "TwistedFermBCParams: Invalid size for boundary. Should be " << Nd << "  but is " << boundary.size() << endl;
      QDP_abort(1);
    }

    if( phases_by_pi.size() != (Nd-1) ) { 
      QDPIO::cerr << "TwistedFermBCParams: Invalid size for phases_by_pi. Should be " << Nd-1 << "  but is " << phases_by_pi.size() << endl;
      QDP_abort(1);
    }

    if( phases_dir.size() != (Nd-1) ) { 
      QDPIO::cerr << "TwistedFermBCParams: Invalid size for phases_dir. Should be " << Nd-1 << "  but is " << phases_dir.size() << endl;
      QDP_abort(1);
    }

    for(int i=0; i < Nd-1; i++) { 
      if( toBool( phases_dir[i] < 0 || phases_dir[i] > Nd-1 ) ) { 
	QDPIO::cerr << "Invalid value in phases_dir, direction " << i << " should be between 0 and " << Nd-1 << " but is " << phases_dir[i] << endl;
      }
    }

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
    write(xml_out, "boundary", param.boundary);
    write(xml_out, "phases_by_pi", param.phases_by_pi);
    write(xml_out, "phases_dir", param.phases_dir);

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
      return new TwistedFermBC<LatticeFermion>(bc.boundary,
					       bc.phases_by_pi,
					       bc.phases_dir);
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
      return new TwistedFermBC< multi1d<LatticeFermion> >(bc.boundary,
					       bc.phases_by_pi,
					       bc.phases_dir);
    }

    //! Name to be used
    const std::string name = TwistedFermBCEnv::name;

    //! Register the fermbc
    const bool registered = TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBC);
  }

}
