// $Id: simple_fermbc_w.cc,v 2.0 2005-09-25 21:04:27 edwards Exp $
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
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      SimpleFermBCParams bc(xml_in, path);
      return new SimpleFermBC<LatticeFermion>(bc.boundary);
    }

    //! Name to be used
    const std::string name = SimpleFermBCEnv::name;

    //! Register the fermbc
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
  }


  //! Name and registration
  namespace WilsonTypeSimpleFermBCArrayEnv
  {
    //! Callback function
    FermBC< multi1d<LatticeFermion> >* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      SimpleFermBCParams bc(xml_in, path);
      return new SimpleFermBC< multi1d<LatticeFermion> >(bc.boundary);
    }

    //! Name to be used
    const std::string name = SimpleFermBCEnv::name;

    //! Register the fermbc
    const bool registered = TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBC);
  }

}
