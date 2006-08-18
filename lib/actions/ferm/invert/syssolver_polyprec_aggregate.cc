// $Id: syssolver_polyprec_aggregate.cc,v 3.2 2006-08-18 15:52:43 edwards Exp $
/*! \file
 *  \brief All PolyPrec system solver constructors
 */

#include "actions/ferm/invert/syssolver_polyprec_aggregate.h"

#include "actions/ferm/invert/syssolver_polyprec_cg.h"

namespace Chroma
{

  //! Registration aggregator
  namespace PolyPrecSysSolverEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PolyPrecSysSolverCGEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
