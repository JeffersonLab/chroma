// $Id: syssolver_mdagm_aggregate.cc,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_cg.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_array.h"

namespace Chroma
{

  //! Registration aggregator
  namespace MdagMSysSolverEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= MdagMSysSolverCGEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  namespace MdagMSysSolverArrayEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= MdagMSysSolverCGArrayEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
