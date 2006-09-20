// $Id: periodic_fermstate_w.cc,v 1.2 2006-09-20 20:31:41 edwards Exp $
/*! \file
 *  \brief Periodic ferm state and a creator
 */

#include "actions/ferm/fermstates/periodic_fermstate_w.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{

  /*! \ingroup fermstates */
  namespace CreatePeriodicFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreatePeriodicFermState<LatticeFermion,
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
	success &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

}

