#include "chromabase.h"
#include "actions/ferm/fermacts/ovext_const_div_by_resp_strategy.h"

using namespace QDP;

namespace Chroma { 
  namespace OvExtConstDivByResPStrategyEnv { 
    
    AbsOvExtTuningStrategy* createStrategy(XMLReader& xml_in,
					   const std::string& path) 
    {
      Real tuning_constant;
      try {
	XMLReader paramtop(xml_in, path);
	read(paramtop, "./TuningConstant", tuning_constant);
      }
      catch(const std::string& e) {
	QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
      
      return new OvExtConstDivByResPStrategy(tuning_constant);
    }
    
    const std::string name = "OVEXT_DIVIDE_BY_CONSTANT_TIMES_RESP";
    const bool registered = TheAbsOvExtTuningStrategyFactory::Instance().registerObject(name, createStrategy);
  }; // end namespace OvExtTuningStrategy
}; // end namespace Chroma
