// -*- C++ -*-
/*! \file
 *  \brief Register MdagM system solvers
 */

#ifndef __multi_syssolver_mdagm_aggregate_h__
#define __multi_syssolver_mdagm_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMMultiSysSolverEnv
  {
    bool registerAll();
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMMultiSysSolverArrayEnv
  {
    bool registerAll();
  }
}

#endif
