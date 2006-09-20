// -*- C++ -*-
// $Id: syssolver_mdagm_aggregate.h,v 3.2 2006-09-20 20:28:00 edwards Exp $
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
