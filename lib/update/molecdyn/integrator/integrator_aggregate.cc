#include "update/molecdyn/integrator/integrator_aggregate.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "update/molecdyn/integrator/lcm_sexton_weingarten.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator.h"

namespace Chroma {

  namespace LCMMDIntegratorAggregateEnv {

    bool registerAll()
    {
      bool success = true; 

      success &= LatColMatPQPLeapfrogIntegratorEnv::registered;
      success &= LatColMatSextonWeingartenIntegratorEnv::registered;
      success &= LatColMatMinimumNorm2IntegratorEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  };

};
