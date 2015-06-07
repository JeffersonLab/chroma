// -*- C++ -*-
/*! \file
 *  \brief Register linop system solvers that solve  PolyPrec*psi=chi
 */

#ifndef __syssolver_polyprec_aggregate_h__
#define __syssolver_polyprec_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup invert */
  namespace PolyPrecSysSolverEnv
  {
    bool registerAll();
  }
}

#endif
