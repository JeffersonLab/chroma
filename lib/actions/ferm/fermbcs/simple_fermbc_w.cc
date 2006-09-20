// $Id: simple_fermbc_w.cc,v 3.1 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/simple_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeSimpleFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml_in, 
						       const std::string& path)
    {
      return new SimpleFermBC<LatticeFermion,
	                      multi1d<LatticeColorMatrix>, 
	                      multi1d<LatticeColorMatrix> >(SimpleFermBCParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "SIMPLE_FERMBC";

    static bool registered = false;

    //! Register all objects
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
