// -*- C++ -*-
// $Id: syssolver_linop.h,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Disambiguator for LinOp system solvers
 */

#ifndef __syssolver_linop_h__
#define __syssolver_linop_h__

#include "linearop.h"
#include "handle.h"
#include "syssolver.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct LinOpSystemSolver : virtual public SystemSolver<T>
  {
  };


  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct LinOpSystemSolverArray : virtual public SystemSolverArray<T>
  {
  };

}


#endif
