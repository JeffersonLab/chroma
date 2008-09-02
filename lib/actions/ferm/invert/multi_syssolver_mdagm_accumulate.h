// -*- C++ -*-
// $Id: multi_syssolver_mdagm_accumulate.h,v 3.1 2008-09-02 20:10:18 bjoo Exp $
/*! \file
 *  \brief Disambiguator for multi-shift MdagM system solvers
 */

#ifndef __multi_syssolver_mdagm_accumulate_h__
#define __multi_syssolver_mdagm_accumulate_h__

#include "linearop.h"
#include "handle.h"
#include "syssolver.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct MdagMMultiSystemSolverAccumulate : virtual public MultiSystemSolverAccumulate<T>
  {
  };


  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct MdagMMultiSystemSolverAccumulateArray : virtual public MultiSystemSolverAccumulateArray<T>
  {
  };

}


#endif
