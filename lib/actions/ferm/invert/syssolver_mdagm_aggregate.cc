// $Id: syssolver_mdagm_aggregate.cc,v 3.2 2006-08-18 15:52:43 edwards Exp $
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
