// $Id: syssolver_polyprec_aggregate.cc,v 3.1 2006-07-03 15:26:09 edwards Exp $
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
