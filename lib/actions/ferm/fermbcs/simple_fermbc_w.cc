// $Id: simple_fermbc_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
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

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
      return foo;
    }

    const bool registered = registerAll();
  }

}
