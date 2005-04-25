#include "chromabase.h"
#include "actions/ferm/fermacts/ovext_neuberger_strategy.h"

using namespace QDP;

namespace Chroma { 
  namespace OvExtNeubergerStrategyEnv { 
    
    AbsOvExtTuningStrategy* createStrategy(XMLReader& xml_in,
					   const std::string& path) 
    {
      return new OvExtNeubergerStrategy();
    }
    
    const std::string name = "OVEXT_NEUBERGER_STRATEGY";
    const bool registered = TheAbsOvExtTuningStrategyFactory::Instance().registerObject(name, createStrategy);
  }; // end namespace OvExtTuningStrategy
}; // end namespace Chroma
