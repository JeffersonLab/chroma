/*! \file
 *  \brief Periodic ferm state and a creator
 */

#include "actions/ferm/fermstates/periodic_fermstate_s.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_s.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_s.h"
#include "actions/ferm/fermbcs/fermbcs_reader_s.h"

namespace Chroma
{
 
  namespace StaggeredCreatePeriodicFermStateEnv 
  { 
    CreateFermState<LatticeStaggeredFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createStag(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreatePeriodicFermState<LatticeStaggeredFermion,
	                                 multi1d<LatticeColorMatrix>, 
	                                 multi1d<LatticeColorMatrix> >();
    }

    const std::string name = "PERIODIC_FERM_STATE";

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

