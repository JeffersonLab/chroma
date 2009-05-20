// $Id: syssolver_linop_aggregate.cc,v 3.13 2009-05-20 15:25:52 bjoo Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "actions/ferm/invert/syssolver_linop_bicgstab.h"
#include "actions/ferm/invert/syssolver_linop_mr.h"
#include "actions/ferm/invert/syssolver_linop_cg_timing.h"
#include "actions/ferm/invert/syssolver_linop_eigcg.h"
#include "actions/ferm/invert/syssolver_linop_richardson_multiprec_clover.h"
#include "actions/ferm/invert/syssolver_linop_rel_bicgstab_clover.h"

#include "actions/ferm/invert/syssolver_linop_cg_array.h"
#include "actions/ferm/invert/syssolver_linop_eigcg_array.h"


namespace Chroma
{

  //! Registration aggregator
  namespace LinOpSysSolverEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// 4D system solvers
	success &= LinOpSysSolverCGEnv::registerAll();
	success &= LinOpSysSolverBiCGStabEnv::registerAll();
	success &= LinOpSysSolverMREnv::registerAll();
	success &= LinOpSysSolverCGTimingEnv::registerAll();
	success &= LinOpSysSolverEigCGEnv::registerAll();
	success &= LinOpSysSolverRichardsonCloverEnv::registerAll();
	success &= LinOpSysSolverReliableBiCGStabCloverEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }


  //! Registration aggregator
  namespace LinOpSysSolverArrayEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// 5D system solvers
	success &= LinOpSysSolverCGArrayEnv::registerAll();
	success &= LinOpSysSolverEigCGArrayEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
