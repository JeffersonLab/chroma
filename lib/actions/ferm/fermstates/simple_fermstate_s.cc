// $Id: simple_fermstate_s.cc,v 1.1 2006-09-19 17:53:37 edwards Exp $
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

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheStaggeredCreateFermStateFactory::Instance().registerObject(name, 
										   createStag);
      return foo;
    }

    const bool registered = registerAll();
  };

}

