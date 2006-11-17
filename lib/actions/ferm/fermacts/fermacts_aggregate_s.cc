// $Id: fermacts_aggregate_s.cc,v 3.1 2006-11-17 02:17:31 edwards Exp $
/*! \file
 *  \brief All Staggered-type fermion actions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_s.h"
#include "actions/ferm/fermacts/fermacts_aggregate_s.h"

#include "actions/ferm/fermacts/asqtad_fermact_s.h"

#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/fermstates/ferm_createstate_aggregate_s.h"


namespace Chroma
{

  //! Registration aggregator
  namespace StaggeredTypeFermActsEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// All system solvers
	success &= LinOpSysSolverEnv::registerAll();
	success &= MdagMSysSolverEnv::registerAll();
	success &= MdagMMultiSysSolverEnv::registerAll();

	// All 4D bcs
	success &= StaggeredTypeFermBCEnv::registerAll();

	// All fermstates
	success &= StaggeredCreateFermStateEnv::registerAll();

	// 4D actions
	success &= AsqtadFermActEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
