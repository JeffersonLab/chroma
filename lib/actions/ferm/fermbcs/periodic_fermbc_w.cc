// $Id: periodic_fermbc_w.cc,v 3.1 2006-09-20 20:28:00 edwards Exp $
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

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
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
