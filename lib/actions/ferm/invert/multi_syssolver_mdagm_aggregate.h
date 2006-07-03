// -*- C++ -*-
// $Id: multi_syssolver_mdagm_aggregate.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
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
    extern const bool registered;
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMMultiSysSolverArrayEnv
  {
    extern const bool registered;
  }
}

#endif
