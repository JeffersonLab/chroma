// $Id: fermbcs_w.cc,v 1.1 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief All fermionic BC
 */

#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include "actions/ferm/fermbcs/simple_fermbc_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeFermBCEnv
  {
    bool registerAll(void) 
    {
      bool success; 
      success = WilsonTypeSimpleFermBCEnv::registered;
      return success;
    }

    const bool registered = registerAll();

    // Helper function for the FermionAction readers
    /*
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeFermion> > reader(XMLReader& xml_in, const std::string& path)
    {
      XMLReader top(xml_in, path);

      bool success = WilsonTypeSimpleFermBCEnv::registered;  // make sure all codes loaded

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
	fermbc = WilsonTypeSimpleFermBCEnv::name;
      }
      else
      {
	QDPIO::cerr << "Error: neither FermionBC group nor boundary found" << endl;
	QDP_abort(1);
      }

      Handle< FermBC<LatticeFermion> > 
	fbc(TheWilsonTypeFermBCFactory::Instance().createObject(fermbc,
								top,
								fermbc_path));

      return fbc;
    }
  }


  //! Name and registration
  namespace WilsonTypeFermBCArrayEnv
  {
    bool registerAll(void) 
    {
      bool success; 
      success = WilsonTypeSimpleFermBCArrayEnv::registered;
      return success;
    }

    const bool registered = registerAll();

    // Helper function for the FermionAction readers
    /*
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC< multi1d<LatticeFermion> > > reader(XMLReader& xml_in, const std::string& path)
    {
      XMLReader top(xml_in, path);

      bool success = WilsonTypeSimpleFermBCEnv::registered;  // make sure all codes loaded

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
	fermbc = WilsonTypeSimpleFermBCArrayEnv::name;
      }
      else
      {
	XMLFileWriter xml_tmp("xml_tmp");
	push(xml_tmp,"XmlOut");
	xml_tmp << top;
	pop(xml_tmp);
	QDPIO::cerr << "path = " << path << endl;
	QDPIO::cerr << "Error: neither FermionBC group nor boundary found" << endl;
	QDP_abort(1);
      }

      Handle< FermBC< multi1d<LatticeFermion> > > 
	fbc(TheWilsonTypeFermBCArrayFactory::Instance().createObject(fermbc,
								     top,
								     fermbc_path));

      return fbc;
    }
  }

}
