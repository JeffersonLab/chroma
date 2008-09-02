// -*- C++ -*-
// $Id: multi_syssolver_mdagm_accumulate_aggregate.h,v 3.1 2008-09-02 20:10:18 bjoo Exp $
/*! \file
 *  \brief Register MdagM system solvers
 */

#ifndef __multi_syssolver_mdagm_accumulate_aggregate_h__
#define __multi_syssolver_mdagm_accumulate_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMMultiSysSolverAccumulateEnv
  {
    bool registerAll();
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMMultiSysSolverAccumulateArrayEnv
  {
    bool registerAll();
  }
}

#endif
