// $Id: multi_syssolver_linop_aggregate.cc,v 1.1 2007-04-11 03:41:36 edwards Exp $
/*! \file
 *  \brief All LinOp system solver constructors
 */

#include "actions/ferm/invert/multi_syssolver_linop_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_linop_mr.h"

namespace Chroma
{
  //! Registration aggregator
  namespace LinOpMultiSysSolverEnv
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
//	success &= LinOpMultiSysSolverMREnv::registerAll();
	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace LinOpMultiSysSolverArrayEnv
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
//	success &= LinOpMultiSysSolverMRArrayEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
