// -*- C++ -*-
/*! \file
 *  \brief Disambiguator for PolyPrec system solvers
 */

#ifndef __syssolver_polyprec_h__
#define __syssolver_polyprec_h__

#include "linearop.h"
#include "handle.h"
#include "syssolver.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct PolyPrecSystemSolver : virtual public SystemSolver<T>
  {
  };

}


#endif
