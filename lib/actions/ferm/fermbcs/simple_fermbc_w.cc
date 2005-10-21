// $Id: simple_fermbc_w.cc,v 2.1 2005-10-21 19:32:16 kostas Exp $
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
#if defined(__APPLE__)
    // fix the gcc bug on MAC OS X
    const std::string name="SIMPLE_FERMBC";
#else
    const std::string name = SimpleFermBCEnv::name;
#endif

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
#if defined(__APPLE__)
    // fix the gcc bug on MAC OS X
    const std::string name="SIMPLE_FERMBC";
#else
    const std::string name = SimpleFermBCEnv::name;
#endif

    //! Register the fermbc
    const bool registered = TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBC);
  }

}
