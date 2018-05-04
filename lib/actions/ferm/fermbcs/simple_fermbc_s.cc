/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/simple_fermbc_s.h"
#include "actions/ferm/fermbcs/fermbc_factory_s.h"

namespace Chroma
{

  //! Name and registration
  namespace StaggeredTypeSimpleFermBCEnv
  {
    //! Callback function
    FermBC<LatticeStaggeredFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      SimpleFermBCParams bc(xml_in, path);
      return new SimpleFermBC<LatticeStaggeredFermion,
	                      multi1d<LatticeColorMatrix>,
	                      multi1d<LatticeColorMatrix> >(bc.boundary);
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
	success &= TheStaggeredTypeFermBCFactory::Instance().registerObject(name, createFermBC);
	registered = true;
      }
      return success;
    }
  }

}
