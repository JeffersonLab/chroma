/*! \file
 *  \brief All MdagM system solver constructors
 */


#include "actions/ferm/invert/syssolver_linop_mrhs_aggregate.h"
#include "actions/ferm/invert/syssolver_mrhs_proxy.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_proxy.h"

namespace Chroma
{

  //! Registration aggregator
  namespace LinOpSysSolverMRHSEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
    	  success &= LinOpSysSolverMRHSProxyEnv::registerAll();
    	  success &= LinOpSysSolverMRHSTwistedProxyEnv::registerAll();
    	  registered = true;
      }
      return success;
    }
  }

}
