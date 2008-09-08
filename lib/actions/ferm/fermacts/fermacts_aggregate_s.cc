// $Id: fermacts_aggregate_s.cc,v 3.4 2008-09-08 20:05:24 bjoo Exp $
/*! \file
 *  \brief All Staggered-type fermion actions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_s.h"
#include "actions/ferm/fermacts/fermacts_aggregate_s.h"

#include "actions/ferm/fermacts/asqtad_fermact_s.h"
#include "actions/ferm/fermacts/hisq_fermact_s.h"
#include "actions/ferm/fermacts/klein_gordon_fermact_s.h"

#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"
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
        success &= MdagMMultiSysSolverAccumulateEnv::registerAll();

	// All 4D bcs
	success &= StaggeredTypeFermBCEnv::registerAll();

	// All fermstates
	success &= StaggeredCreateFermStateEnv::registerAll();

	// 4D actions
	success &= AsqtadFermActEnv::registerAll();
	success &= HisqFermActEnv::registerAll();
	success &= KleinGordonFermActEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
