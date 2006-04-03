// $Id: periodic_fermbc_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
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
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml_in, 
						       const std::string& path)
    {
      return new PeriodicFermBC<LatticeFermion, 
	                        multi1d<LatticeColorMatrix>, 
                                multi1d<LatticeColorMatrix> >();
    }

    //! Name to be used
    const std::string name = "PERIODIC_FERMBC";

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
