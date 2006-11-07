// -*- C++ -*-
// $Id: integrator_aggregate.h,v 3.3 2006-11-07 23:11:01 bjoo Exp $

#ifndef INTEGRATOR_AGGREGATE_H
#define INTEGRATOR_AGGREGATE_H


#include "chromabase.h"

namespace Chroma 
{ 
  //! Registration aggregator
  /*! @ingroup monomial */
  namespace LCMMDIntegratorAggregateEnv 
  {
    bool registerAll();
  }

  namespace LCMMDComponentIntegratorAggregateEnv 
  {
    bool registerAll();
  }
}


#endif
