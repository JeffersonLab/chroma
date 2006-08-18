// $Id: syssolver_linop_aggregate.cc,v 3.2 2006-08-18 15:52:43 edwards Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "actions/ferm/invert/syssolver_linop_cg_array.h"

namespace Chroma
{

  //! Registration aggregator
  namespace LinOpSysSolverEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // 4D system solvers
      success &= LinOpSysSolverCGEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  namespace LinOpSysSolverArrayEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // 5D system solvers
      success &= LinOpSysSolverCGArrayEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
