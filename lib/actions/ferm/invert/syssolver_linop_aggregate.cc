// $Id: syssolver_linop_aggregate.cc,v 3.4 2006-11-16 20:39:15 bjoo Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "actions/ferm/invert/syssolver_linop_cg_timing.h"
#include "actions/ferm/invert/syssolver_linop_cg_array.h"

namespace Chroma
{

  //! Registration aggregator
  namespace LinOpSysSolverEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// 4D system solvers
	success &= LinOpSysSolverCGEnv::registerAll();
	success &= LinOpSysSolverCGTimingEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace LinOpSysSolverArrayEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// 5D system solvers
	success &= LinOpSysSolverCGArrayEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
