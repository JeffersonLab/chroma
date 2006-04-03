// -*- C++ -*-
// $Id: ovext_tuning_strategy.h,v 3.0 2006-04-03 04:58:45 edwards Exp $
/*! \file
 *  \brief Ovext tuning strategy
 */

#ifndef ovext_tuning_strategy_h
#define ovext_tuning_strategy_h

#include "singleton.h"
#include "objfactory.h"

#include "chromabase.h"

using namespace QDP;

namespace Chroma 
{ 
  // Base class for tuning strategies. Only one public function: operator()

  //! Ovext tuning strategy
  /*! @ingroup fermacts */
  class AbsOvExtTuningStrategy 
  { 
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


  /*! @ingroup fermacts */
  typedef SingletonHolder< 
    ObjectFactory< AbsOvExtTuningStrategy,
		   std::string,
		   TYPELIST_2(XMLReader&, const std::string&),
		   AbsOvExtTuningStrategy* (*)(XMLReader&, const std::string&),
		   StringFactoryError> > 
  TheAbsOvExtTuningStrategyFactory;
  
}


#endif
