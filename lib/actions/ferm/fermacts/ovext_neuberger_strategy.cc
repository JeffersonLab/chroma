// -*- C++ -*-
// $Id: ovext_neuberger_strategy.cc,v 3.0 2006-04-03 04:58:45 edwards Exp $
/*! \file
 *  \brief Ovext Neuberger rescale strategy
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/ovext_neuberger_strategy.h"

using namespace QDP;

namespace Chroma 
{ 

  namespace OvExtNeubergerStrategyEnv 
  { 
    AbsOvExtTuningStrategy* createStrategy(XMLReader& xml_in,
					   const std::string& path) 
    {
      return new OvExtNeubergerStrategy();
    }
    
    const std::string name = "OVEXT_NEUBERGER_STRATEGY";
    const bool registered = TheAbsOvExtTuningStrategyFactory::Instance().registerObject(name, createStrategy);
  } // end namespace OvExtTuningStrategy

} // end namespace Chroma
