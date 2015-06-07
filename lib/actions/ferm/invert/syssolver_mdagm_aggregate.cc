/*! \file
 *  \brief All MdagM system solver constructors
 */

#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"


#include "actions/ferm/invert/syssolver_mdagm_cg.h"
#include "actions/ferm/invert/syssolver_mdagm_bicgstab.h"
#include "actions/ferm/invert/syssolver_mdagm_ibicgstab.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_timing.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_array.h"
#include "actions/ferm/invert/syssolver_mdagm_eigcg.h"
#include "actions/ferm/invert/syssolver_mdagm_richardson_multiprec_clover.h"
#include "actions/ferm/invert/syssolver_mdagm_rel_bicgstab_clover.h"
#include "actions/ferm/invert/syssolver_mdagm_rel_ibicgstab_clover.h"
#include "actions/ferm/invert/syssolver_mdagm_rel_cg_clover.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_lf_clover.h"

#include "chroma_config.h"
#ifdef BUILD_QUDA
#include "actions/ferm/invert/quda_solvers/syssolver_mdagm_clover_quda_w.h"
#include "actions/ferm/invert/quda_solvers/syssolver_mdagm_wilson_quda_w.h"
#endif

#ifdef BUILD_QPHIX
#include "actions/ferm/invert/qphix/syssolver_mdagm_clover_qphix_w.h"
#endif

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
	success &= MdagMSysSolverIBiCGStabEnv::registerAll();
	success &= MdagMSysSolverEigCGEnv::registerAll();
	success &= MdagMSysSolverRichardsonCloverEnv::registerAll();
	success &= MdagMSysSolverReliableBiCGStabCloverEnv::registerAll();
	success &= MdagMSysSolverReliableIBiCGStabCloverEnv::registerAll();
	success &= MdagMSysSolverReliableCGCloverEnv::registerAll();
	success &= MdagMSysSolverCGLFCloverEnv::registerAll();
#ifdef BUILD_QUDA
	success &= MdagMSysSolverQUDACloverEnv::registerAll();
	success &= MdagMSysSolverQUDAWilsonEnv::registerAll();
#endif

#ifdef BUILD_QPHIX
	success &= MdagMSysSolverQPhiXCloverEnv::registerAll();
#endif
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
