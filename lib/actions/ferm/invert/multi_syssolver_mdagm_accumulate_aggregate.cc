// $Id: multi_syssolver_mdagm_accumulate_aggregate.cc,v 3.2 2008-09-06 18:35:35 bjoo Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg_accumulate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_cg_accumulate_array.h"

namespace Chroma
{
  //! Registration aggregator
  namespace MdagMMultiSysSolverAccumulateEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Sources
	success &= MdagMMultiSysSolverAccumulateCGEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

  /* This doesnt work yet */

  //! Registration aggregator
  namespace MdagMMultiSysSolverAccumulateArrayEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Sources
	success &= MdagMMultiSysSolverCGAccumulateArrayEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }


}
