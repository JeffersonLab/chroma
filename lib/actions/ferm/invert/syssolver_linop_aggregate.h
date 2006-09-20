// -*- C++ -*-
// $Id: syssolver_linop_aggregate.h,v 3.2 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Register linop system solvers that solve  M*psi=chi
 */

#ifndef __syssolver_linop_aggregate_h__
#define __syssolver_linop_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup invert */
  namespace LinOpSysSolverEnv
  {
    bool registerAll();
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace LinOpSysSolverArrayEnv
  {
    bool registerAll();
  }
}

#endif
