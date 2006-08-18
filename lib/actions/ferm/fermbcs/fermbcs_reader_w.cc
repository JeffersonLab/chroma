// $Id: fermbcs_reader_w.cc,v 3.1 2006-08-18 15:52:43 edwards Exp $
/*! \file
 *  \brief Fermionic BC reader
 */

#include "actions/ferm/fermbcs/fermbcs_reader_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include "actions/ferm/fermbcs/simple_fermbc_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeFermBCEnv
  {
    // Helper function for the FermionAction readers
    /*
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeFermion,
		   multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
							  const std::string& path)
    {
      XMLReader top(xml_in, path);

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

      Handle< FermBC<LatticeFermion,
	             multi1d<LatticeColorMatrix>,
	             multi1d<LatticeColorMatrix> > > 
	fbc(TheWilsonTypeFermBCFactory::Instance().createObject(fermbc,
								top,
								fermbc_path));

      return fbc;
    }
  }

}
