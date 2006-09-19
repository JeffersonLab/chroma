// $Id: ferm_createstate_reader_w.cc,v 1.1 2006-09-19 17:53:37 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#include "chromabase.h"

#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/simple_fermstate_w.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"

namespace Chroma
{

  // State reader
  namespace CreateFermStateEnv
  {

    // Helper function for the FermAction readers
    Handle< CreateFermState<LatticeFermion,
			    multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								   const std::string& path)
    {
      XMLReader top(xml_in, path);

      std::string fermstate;
      std::string fermstate_path;
      if (top.count("FermState") != 0)
      {
	fermstate_path = "FermState";
	read(top, fermstate_path + "/Name", fermstate);
      }
      else
      {
//	QDPIO::cerr << "Error: FermState group not found" << endl;
//	QDP_abort(1);

	fermstate_path = ".";
	fermstate = Chroma::CreateSimpleFermStateEnv::name;
      }

      Handle< CreateFermState<
	LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	cgs(TheCreateFermStateFactory::Instance().createObject(fermstate,
							       top,
							       fermstate_path));

      return cgs;
    }

  }

} // end chroma namespace
