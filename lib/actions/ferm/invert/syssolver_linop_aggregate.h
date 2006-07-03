// -*- C++ -*-
// $Id: syssolver_linop_aggregate.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
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
    extern const bool registered;
  }


  //! Registration aggregator
  /*! @ingroup invert */
  namespace LinOpSysSolverArrayEnv
  {
    extern const bool registered;
  }
}

#endif
