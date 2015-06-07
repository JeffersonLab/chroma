// -*- C++ -*-
/*! \file
 *  \brief Register MdagM system solvers
 */

#ifndef __syssolver_mdagm_aggregate_h__
#define __syssolver_mdagm_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMSysSolverEnv
  {
    bool registerAll();
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMSysSolverArrayEnv
  {
    bool registerAll();
  }
}

#endif
