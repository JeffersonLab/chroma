#include "update/molecdyn/integrator/integrator_aggregate.h"
#include "update/molecdyn/integrator/lcm_exp_sdt.h"
#include "update/molecdyn/integrator/lcm_exp_tdt.h"
#include "update/molecdyn/integrator/lcm_sts_leapfrog_recursive.h"
#include "update/molecdyn/integrator/lcm_tst_leapfrog_recursive.h"
#include "update/molecdyn/integrator/lcm_sts_min_norm2_recursive.h"
#include "update/molecdyn/integrator/lcm_sts_min_norm2_recursive_dtau.h"
#include "update/molecdyn/integrator/lcm_tst_min_norm2_recursive.h"
#include "update/molecdyn/integrator/lcm_tst_min_norm2_recursive_dtau.h"
#include "update/molecdyn/integrator/lcm_4mn5fv_recursive.h"
#include "update/molecdyn/integrator/lcm_4mn5fp_recursive.h"
#include "update/molecdyn/integrator/lcm_4mn4fp_recursive.h"
#include "update/molecdyn/integrator/lcm_creutz_gocksch_4_recursive.h"
namespace Chroma 
{

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
	success &=  LatColMatTSTLeapfrogRecursiveIntegratorEnv::registerAll();
	success &=  LatColMatSTSMinNorm2RecursiveIntegratorEnv::registerAll();
	success &=  LatColMatSTSMinNorm2DTauRecursiveIntegratorEnv::registerAll();
	success &=  LatColMatTSTMinNorm2RecursiveIntegratorEnv::registerAll();
	success &=  LatColMatTSTMinNorm2DTauRecursiveIntegratorEnv::registerAll();

	success &=  LatColMat4MN4FPRecursiveIntegratorEnv::registerAll();
	success &=  LatColMat4MN5FVRecursiveIntegratorEnv::registerAll();
	success &=  LatColMat4MN5FPRecursiveIntegratorEnv::registerAll();
	success &=  LatColMatCreutzGocksch4RecursiveIntegratorEnv::registerAll();
	registered = true;
      }
      return success;
    }

  }
}
