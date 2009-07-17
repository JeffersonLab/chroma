// -*- C++ -*-
// $Id: syssolver_linop.h,v 3.2 2009-07-17 19:14:46 bjoo Exp $
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
  
  /* NB: Previously these were declared as 'virtual public SystemSolver<T>'
     BUT That seemed to break XLC in a Bad Way */
  template<typename T>
  class LinOpSystemSolver : public SystemSolver<T>
  {
  };


  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  
  /* NB: Previously this was declared as virtual public SystemSolverArray<T> 
     but that broke the xlC build */
  template<typename T>
  class LinOpSystemSolverArray : public SystemSolverArray<T>
  {
  };

}


#endif
