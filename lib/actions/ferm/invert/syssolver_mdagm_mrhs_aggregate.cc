/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

namespace Chroma
{

  //! Registration aggregator
  namespace MdagMSysSolverMRHSEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
    	  registered = true;
      }
      return success;
    }
  }
}
