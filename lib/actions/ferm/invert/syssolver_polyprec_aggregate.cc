// $Id: syssolver_polyprec_aggregate.cc,v 3.3 2006-09-20 20:28:00 edwards Exp $
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
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Sources
	success &= PolyPrecSysSolverCGEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
