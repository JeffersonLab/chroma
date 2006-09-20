// -*- C++ -*-
// $Id: ovext_neuberger_strategy.cc,v 3.1 2006-09-20 20:27:59 edwards Exp $
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

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheAbsOvExtTuningStrategyFactory::Instance().registerObject(name, createStrategy);
	registered = true;
      }
      return success;
    }

  } // end namespace OvExtTuningStrategy
} // end namespace Chroma
