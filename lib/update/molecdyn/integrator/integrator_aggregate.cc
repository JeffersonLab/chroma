#include "update/molecdyn/integrator/integrator_aggregate.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog.h"
#include "update/molecdyn/integrator/lcm_sexton_weingarten.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator.h"
#include "update/molecdyn/integrator/lcm_sw_min_mixed.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_integrator_mts.h"
#include "update/molecdyn/integrator/lcm_minimum_norm2_qpq_integrator_mts.h"
#include "update/molecdyn/integrator/lcm_pqp_leapfrog_mts.h"

#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "update/molecdyn/integrator/lcm_exp_tdt.h"
#include "update/molecdyn/integrator/lcm_sts_leapfrog_recursive.h"
#include "update/molecdyn/integrator/lcm_min_norm2_recursive.h"
#include "update/molecdyn/integrator/lcm_4mn5fv_recursive.h"
#include "update/molecdyn/integrator/lcm_4mn4fp_recursive.h"
#include "update/molecdyn/integrator/lcm_creutz_gocksch_4_recursive.h"
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
	//	success &= LatColMatSexWeinMixedIntegratorEnv::registerAll();
	success &= LatColMatMinimumNorm2IntegratorMtsEnv::registerAll();
	success &= LatColMatMinimumNorm2QPQIntegratorMtsEnv::registerAll();
	success &= LatColMatPQPLeapfrogIntegratorMtsEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

  namespace LCMMDComponentIntegratorAggregateEnv 
  {
    namespace 
    {
      //! Local registration flag
      bool registered = false;
    }

    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &=  LatColMatExpSdtIntegratorEnv::registerAll();
	success &=  LatColMatExpTdtIntegratorEnv::registerAll();
	success &=  LatColMatSTSLeapfrogRecursiveIntegratorEnv::registerAll();
	success &=  LatColMatMinNorm2RecursiveIntegratorEnv::registerAll();

	success &=  LatColMat4MN4FPRecursiveIntegratorEnv::registerAll();
	success &=  LatColMat4MN5FVRecursiveIntegratorEnv::registerAll();
	success &=  LatColMatCreutzGocksch4RecursiveIntegratorEnv::registerAll();
	registered = true;
      }
      return success;
    }

  }
}
