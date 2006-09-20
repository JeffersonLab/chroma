#include "update/molecdyn/integrator/integrator_aggregate.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "update/molecdyn/integrator/lcm_sexton_weingarten.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator.h"
#include "update/molecdyn/integrator/lcm_sw_min_mixed.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator_mts.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_qpq_integrator_mts.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog_mts.h"

namespace Chroma 
{

  namespace LCMMDIntegratorAggregateEnv 
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= LatColMatPQPLeapfrogIntegratorEnv::registerAll();
	success &= LatColMatSextonWeingartenIntegratorEnv::registerAll();
	success &= LatColMatMinimumNorm2IntegratorEnv::registerAll();
	success &= LatColMatSexWeinMixedIntegratorEnv::registerAll();
	success &= LatColMatMinimumNorm2IntegratorMtsEnv::registerAll();
	success &= LatColMatMinimumNorm2QPQIntegratorMtsEnv::registerAll();
	success &= LatColMatPQPLeapfrogIntegratorMtsEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
