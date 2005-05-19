#include "update/molecdyn/integrator/integrator_aggregate.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "update/molecdyn/integrator/lcm_sexton_weingarten.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator.h"
#include "update/molecdyn/integrator/lcm_sw_min_mixed.h"
namespace Chroma {

  namespace LCMMDIntegratorAggregateEnv {

    bool registerAll()
    {
      bool success = true; 

      success &= LatColMatPQPLeapfrogIntegratorEnv::registered;
      success &= LatColMatSextonWeingartenIntegratorEnv::registered;
      success &= LatColMatMinimumNorm2IntegratorEnv::registered;
      success &= LatColMatSexWeinMixedIntegratorEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  };

};
