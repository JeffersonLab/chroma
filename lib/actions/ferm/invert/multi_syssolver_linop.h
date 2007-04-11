// -*- C++ -*-
// $Id: multi_syssolver_linop.h,v 1.1 2007-04-11 03:41:36 edwards Exp $
/*! \file
 *  \brief Disambiguator for multi-shift linop system solvers
 */

#ifndef __multi_syssolver_linop_h__
#define __multi_syssolver_linop_h__

#include "linearop.h"
#include "handle.h"
#include "syssolver.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct LinOpMultiSystemSolver : virtual public MultiSystemSolver<T>
  {
  };


  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct LinOpMultiSystemSolverArray : virtual public MultiSystemSolverArray<T>
  {
  };

}


#endif
