/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/syssolver_mrhs_proxy.h"
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
    	  MdagMSysSolverMRHSProxyEnv::registerAll();
    	  registered = true;
      }
      return success;
    }
  }
}
