// $Id: predictor_aggregate.cc,v 1.2 2005-02-23 22:24:15 bjoo Exp $
/*! \file
 *  \brief Chrono predictor aggregator
 */

#include "update/molecdyn/predictor/predictor_aggregate.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"
#include "update/molecdyn/predictor/last_solution_predictor.h"
#include "update/molecdyn/predictor/linear_extrap_predictor.h"

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
