// $Id: periodic_fermbc_w.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Periodic fermionic BC
 */

#include "actions/ferm/fermbcs/periodic_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypePeriodicFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      return new PeriodicFermBC<LatticeFermion>();
    }

    //! Callback function
    FermBC< multi1d<LatticeFermion> >* createFermBCArray(XMLReader& xml_in, const std::string& path)
    {
      return new PeriodicFermBC< multi1d<LatticeFermion> >();
    }

    //! Name to be used
    const std::string name = "PERIODIC_FERMBC";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
      foo &= Chroma::TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBCArray);
      return foo;
    }

    const bool registered = registerAll();
  }

}
