// $Id: syssolver_linop_aggregate.cc,v 3.21 2009-11-02 21:34:08 kostas Exp $
/*! \file
 *  \brief All MdagM system solver constructors
 */


#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "actions/ferm/invert/syssolver_linop_bicgstab.h"
#include "actions/ferm/invert/syssolver_linop_ibicgstab.h"
#include "actions/ferm/invert/syssolver_linop_bicrstab.h"
#include "actions/ferm/invert/syssolver_linop_mr.h"
#include "actions/ferm/invert/syssolver_linop_cg_timing.h"
#include "actions/ferm/invert/syssolver_linop_eigcg.h
#include "actions/ferm/invert/syssolver_linop_OPTeigbicg.h"
#include "actions/ferm/invert/syssolver_linop_richardson_multiprec_clover.h"
#include "actions/ferm/invert/syssolver_linop_rel_bicgstab_clover.h"
#include "actions/ferm/invert/syssolver_linop_rel_ibicgstab_clover.h"
#include "actions/ferm/invert/syssolver_linop_rel_cg_clover.h"


#include "chroma_config.h"
#ifdef BUILD_QUDA
#include "actions/ferm/invert/quda_solvers/syssolver_linop_quda_wilson.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_quda_clover.h"
#endif

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
	success &= LinOpSysSolverBiCRStabEnv::registerAll();
	success &= LinOpSysSolverIBiCGStabEnv::registerAll();
	success &= LinOpSysSolverMREnv::registerAll();
	success &= LinOpSysSolverCGTimingEnv::registerAll();
	success &= LinOpSysSolverEigCGEnv::registerAll();
	success &= LinOpSysSolverRichardsonCloverEnv::registerAll();
	success &= LinOpSysSolverReliableBiCGStabCloverEnv::registerAll();
	success &= LinOpSysSolverReliableIBiCGStabCloverEnv::registerAll();
	success &= LinOpSysSolverReliableCGCloverEnv::registerAll();
     
#ifdef BUILD_QUDA
#warning Registering QUDA
	success &=LinOpSysSolverQUDAWilsonEnv::registerAll();
	success &=LinOpSysSolverQUDACloverEnv::registerAll();
#endif
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
