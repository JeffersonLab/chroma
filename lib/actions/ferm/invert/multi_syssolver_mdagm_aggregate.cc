// $Id: multi_syssolver_mdagm_aggregate.cc,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_cg_array.h"

namespace Chroma
{
  //! Registration aggregator
  namespace MdagMMultiSysSolverEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= MdagMMultiSysSolverCGEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  namespace MdagMMultiSysSolverArrayEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= MdagMMultiSysSolverCGArrayEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
