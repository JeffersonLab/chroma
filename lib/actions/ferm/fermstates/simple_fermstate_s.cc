// $Id: simple_fermstate_s.cc,v 1.2 2006-09-20 20:31:41 edwards Exp $
/*! \file
 *  \brief Simple ferm state and a creator
 */

#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/simple_fermstate_s.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_s.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_s.h"
#include "actions/ferm/fermbcs/fermbcs_reader_s.h"

namespace Chroma
{
 
  namespace StaggeredCreateSimpleFermStateEnv 
  { 
    CreateFermState<LatticeStaggeredFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createStag(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreateSimpleFermState<LatticeStaggeredFermion,
	                               multi1d<LatticeColorMatrix>, 
	                               multi1d<LatticeColorMatrix> >(StaggeredTypeFermBCEnv::reader(xml, 
												    path));
    }

    const std::string name = "SIMPLE_FERM_STATE";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheStaggeredCreateFermStateFactory::Instance().registerObject(name, createStag);
	registered = true;
      }
      return success;
    }
  }

}

