// -*- C++ -*-
/*! \file
 *  \brief Disambiguator for multi-shift MdagM system solvers
 */

#ifndef __multi_syssolver_mdagm_h__
#define __multi_syssolver_mdagm_h__

#include "linearop.h"
#include "handle.h"
#include "syssolver.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct MdagMMultiSystemSolver : virtual public MultiSystemSolver<T>
  {
  };

  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct MdagMMultiSystemSolverArray : virtual public MultiSystemSolverArray<T>
  {
  };

}


#endif
