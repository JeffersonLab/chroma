// $Id: multi_syssolver_mdagm_aggregate.cc,v 3.4 2008-09-06 18:35:35 bjoo Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_cg_array.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_cg_chrono_clover.h"

#include "chroma_config.h"
#ifdef BUILD_QUDA
#include "actions/ferm/invert/quda_solvers/multi_syssolver_mdagm_cg_clover_quda_w.h"
#include "actions/ferm/invert/quda_solvers/multi_syssolver_mdagm_cg_wilson_quda_w.h"
#endif

namespace Chroma
{
  //! Registration aggregator
  namespace MdagMMultiSysSolverEnv
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
	success &= MdagMMultiSysSolverCGEnv::registerAll();
	success &= MdagMMultiSysSolverCGChronoCloverEnv::registerAll();
#ifdef BUILD_QUDA
	success &= MdagMMultiSysSolverCGQudaCloverEnv::registerAll();
	success &= MdagMMultiSysSolverCGQudaWilsonEnv::registerAll();
#endif

	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace MdagMMultiSysSolverArrayEnv
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
	success &= MdagMMultiSysSolverCGArrayEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
