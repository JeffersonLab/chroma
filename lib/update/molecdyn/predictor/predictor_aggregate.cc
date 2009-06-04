// $Id: predictor_aggregate.cc,v 3.4 2009-06-04 20:29:13 bjoo Exp $
/*! \file
 *  \brief Chrono predictor aggregator
 */

#include "update/molecdyn/predictor/predictor_aggregate.h"

#include "update/molecdyn/predictor/null_predictor.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/predictor/last_solution_predictor.h"
#include "update/molecdyn/predictor/linear_extrap_predictor.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"
#include "update/molecdyn/predictor/mre_initcg_extrap_predictor.h"

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
	success &= LastSolution5DChronoPredictorEnv::registerAll();
	success &= LinearExtrapolation4DChronoPredictorEnv::registerAll();
	success &= LinearExtrapolation5DChronoPredictorEnv::registerAll();
	success &= MinimalResidualExtrapolation4DChronoPredictorEnv::registerAll();
	success &= MinimalResidualExtrapolation5DChronoPredictorEnv::registerAll();
	success &= MREInitCG4DChronoPredictorEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
