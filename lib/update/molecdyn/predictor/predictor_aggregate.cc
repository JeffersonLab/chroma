// $Id: predictor_aggregate.cc,v 1.3 2005-02-23 23:48:10 bjoo Exp $
/*! \file
 *  \brief Chrono predictor aggregator
 */

#include "update/molecdyn/predictor/predictor_aggregate.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/predictor/last_solution_predictor.h"
#include "update/molecdyn/predictor/linear_extrap_predictor.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"

namespace Chroma
{

  //! Name and registration
  namespace ChronoPredictorAggregrateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= ZeroGuess4DChronoPredictorEnv::registered;
      success &= ZeroGuess5DChronoPredictorEnv::registered;
      success &= LastSolution4DChronoPredictorEnv::registered;  
      success &= LastSolution5DChronoPredictorEnv::registered;
      success &= LinearExtrapolation4DChronoPredictorEnv::registered;
      success &= LinearExtrapolation5DChronoPredictorEnv::registered;
      success &= MinimalResidualExtrapolation4DChronoPredictorEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
