// $Id: fermbcs_s.cc,v 1.2 2005-01-12 04:44:19 edwards Exp $
/*! \file
 *  \brief All fermionic BC
 */

#include "actions/ferm/fermbcs/fermbc_factory_s.h"
#include "actions/ferm/fermbcs/simple_fermbc_s.h"

namespace Chroma
{

  //! Name and registration
  namespace StaggeredTypeFermBCEnv
  {
    bool registerAll(void) 
    {
      bool success; 
      success = StaggeredTypeSimpleFermBCEnv::registered;
      return success;
    }

    const bool registered = registerAll();


    // Helper function for the FermionAction readers
    /*
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeStaggeredFermion> > reader(XMLReader& xml_in, const std::string& path)
    {
      XMLReader top(xml_in, path);

      bool success = registered;  // make sure all codes loaded

      std::string fermbc;
      std::string fermbc_path;
      if (top.count("FermionBC") != 0)
      {
	fermbc_path = "FermionBC";
	read(top, fermbc_path + "/FermBC", fermbc);
      }
      else if (top.count("boundary") != 0)
      {
	fermbc_path = ".";
	fermbc = StaggeredTypeSimpleFermBCEnv::name;
      }
      else
      {
	QDPIO::cerr << "Error: neither FermionBC group nor boundary found" << endl;
	QDP_abort(1);
      }

      Handle< FermBC<LatticeStaggeredFermion> > 
	fbc(TheStaggeredTypeFermBCFactory::Instance().createObject(fermbc,
								   top,
								   fermbc_path));

      return fbc;
    }
  }

}
