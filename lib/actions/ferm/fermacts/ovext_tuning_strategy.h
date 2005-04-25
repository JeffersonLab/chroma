#ifndef ovext_tuning_strategy_h
#define ovext_tuning_strategy_h

#include "singleton.h"
#include "objfactory.h"

#include "chromabase.h"

using namespace QDP;

namespace Chroma { 
    // Base class for tuning strategies. Only one public function: operator()

    class AbsOvExtTuningStrategy { 
    public:
      // virtual destructor
      virtual ~AbsOvExtTuningStrategy(void) {}
  
      // operator()
      virtual void operator()(multi1d<Real>& beta,
		      const Real& coeffP,
		      const multi1d<Real>& resP,
		      const multi1d<Real>& resQ,
		      const Real& Mass) const = 0;
    };




    typedef SingletonHolder< 
      ObjectFactory< AbsOvExtTuningStrategy,
      std::string,
      TYPELIST_2(XMLReader&, const std::string&),
      AbsOvExtTuningStrategy* (*)(XMLReader&, const std::string&),
      StringFactoryError> > TheAbsOvExtTuningStrategyFactory;

};


#endif
