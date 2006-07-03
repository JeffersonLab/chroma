// -*- C++ -*-
// $Id: syssolver_polyprec_aggregate.h,v 3.1 2006-07-03 15:26:09 edwards Exp $
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
    extern const bool registered;
  }
}

#endif
