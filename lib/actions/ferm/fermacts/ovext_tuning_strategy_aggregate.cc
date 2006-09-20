// -*- C++ -*-
// $Id: ovext_tuning_strategy_aggregate.cc,v 3.1 2006-09-20 20:27:59 edwards Exp $
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
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= OvExtConstantStrategyEnv::registerAll();
	success &= OvExtConstDivByResPStrategyEnv::registerAll();
	success &= OvExtNeubergerStrategyEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }
}
