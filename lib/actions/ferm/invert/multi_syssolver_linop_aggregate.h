// -*- C++ -*-
// $Id: multi_syssolver_linop_aggregate.h,v 1.1 2007-04-11 03:41:36 edwards Exp $
/*! \file
 *  \brief Register LinOp system solvers
 */

#ifndef __multi_syssolver_linop_aggregate_h__
#define __multi_syssolver_linop_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup invert */
  namespace LinOpMultiSysSolverEnv
  {
    bool registerAll();
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace LinOpMultiSysSolverArrayEnv
  {
    bool registerAll();
  }
}

#endif
