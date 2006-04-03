// -*- C++ -*-
// $Id: ovext_tuning_strategy_aggregate.cc,v 3.0 2006-04-03 04:58:45 edwards Exp $
/*! \file
 *  \brief Ovext tuning strategy
 */

#include "ovext_tuning_strategy_aggregate.h"
#include "actions/ferm/fermacts/ovext_constant_strategy.h"
#include "actions/ferm/fermacts/ovext_const_div_by_resp_strategy.h"
#include "actions/ferm/fermacts/ovext_neuberger_strategy.h"

namespace Chroma 
{ 
  namespace OvExtTuningStrategyAggregateEnv 
  {
    bool registerAll()
    {
      bool success = true;
      success &= OvExtConstantStrategyEnv::registered;
      success &= OvExtConstDivByResPStrategyEnv::registered;
      success &= OvExtNeubergerStrategyEnv::registered;
      return success;
    }
    
    const bool registered = registerAll();
  };
}
