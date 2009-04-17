// $Id: syssolver_mdagm_aggregate.cc,v 3.7 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"


#include "actions/ferm/invert/syssolver_mdagm_cg.h"
#include "actions/ferm/invert/syssolver_mdagm_bicgstab.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_timing.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_array.h"
#include "actions/ferm/invert/syssolver_mdagm_eigcg.h"
#include "actions/ferm/invert/syssolver_mdagm_richardson_multiprec_clover.h"

namespace Chroma
{

  //! Registration aggregator
  namespace MdagMSysSolverEnv
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
	success &= MdagMSysSolverCGEnv::registerAll();
	success &= MdagMSysSolverCGTimingsEnv::registerAll();
	success &= MdagMSysSolverBiCGStabEnv::registerAll();
	success &= MdagMSysSolverEigCGEnv::registerAll();
	success &= MdagMSysSolverRichardsonCloverEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace MdagMSysSolverArrayEnv
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
	success &= MdagMSysSolverCGArrayEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
