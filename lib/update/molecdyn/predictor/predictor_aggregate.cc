/*! \file
 *  \brief Chrono predictor aggregator
 */

#include "chroma_config.h"
#include "update/molecdyn/predictor/predictor_aggregate.h"

#include "update/molecdyn/predictor/null_predictor.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/predictor/last_solution_predictor.h"
//Include MG predictor.
#include "update/molecdyn/predictor/MG_predictor.h"
#include "update/molecdyn/predictor/linear_extrap_predictor.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"
#include "update/molecdyn/predictor/mre_initcg_extrap_predictor.h"

#ifdef BUILD_QUDA
#include "update/molecdyn/predictor/quda_predictor.h"
#endif

namespace Chroma
{

  //! Name and registration
  namespace ChronoPredictorAggregrateEnv
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
	success &= Null4DChronoPredictorEnv::registerAll();
	success &= Null5DChronoPredictorEnv::registerAll();
	success &= ZeroGuess4DChronoPredictorEnv::registerAll();
	success &= ZeroGuess5DChronoPredictorEnv::registerAll();
	success &= LastSolution4DChronoPredictorEnv::registerAll();
	//Added special MG predictor!
	success &= MG4DChronoPredictorEnv::registerAll();  
	success &= LastSolution5DChronoPredictorEnv::registerAll();
	success &= LinearExtrapolation4DChronoPredictorEnv::registerAll();
	success &= LinearExtrapolation5DChronoPredictorEnv::registerAll();
	success &= MinimalResidualExtrapolation4DChronoPredictorEnv::registerAll();
	success &= MinimalResidualExtrapolation5DChronoPredictorEnv::registerAll();

#if ! defined (QDP_IS_QDPJIT2)
	success &= MREInitCG4DChronoPredictorEnv::registerAll();
#endif

#ifdef BUILD_QUDA
	success &= QUDA4DChronoPredictorEnv::registerAll();
#endif
	registered = true;
      }
      return success;
    }
  }

}
