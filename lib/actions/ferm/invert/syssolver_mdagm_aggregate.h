// -*- C++ -*-
// $Id: syssolver_mdagm_aggregate.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
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
    extern const bool registered;
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace MdagMSysSolverArrayEnv
  {
    extern const bool registered;
  }
}

#endif
