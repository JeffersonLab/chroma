// $Id: schr1link_fermbc_w.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 * @brief 1-Link Schroedinger function boundary conditions
 */

#if 0


#include "actions/ferm/fermbcs/schr1link_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include "actions/ferm/fermbcs/schr_fermbc_params.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeSchr1LinkFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      Schr1LinkFermBCParams bc(xml_in, path);      
      return new Schr1LinkFermBC<LatticeFermion>(bc.theta);
    }

    //! Name to be used
    const std::string name = "SCHROEDINGER_1LINK_FERMBC";

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
