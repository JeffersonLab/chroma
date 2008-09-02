// $Id: multi_syssolver_mdagm_accumulate_aggregate.cc,v 3.1 2008-09-02 20:10:18 bjoo Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg_accumulate.h"

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
#if 0
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
#endif

}
