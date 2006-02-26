// $Id: schr2link_fermbc_w.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 * @brief 2-Link Schroedinger function boundary conditions
 */

#if 0

#include "actions/ferm/fermbcs/schr2link_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeSchr2linkFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      return new Schr2linkFermBC<LatticeFermion>(Schr2linkFermBCParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "SCHROEDINGER_2LINK_FERMBC";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
//      foo &= Chroma::TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBCArray);
      return foo;
    }

    const bool registered = registerAll();
  }

}

#endif
