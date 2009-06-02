// -*- C++ -*-
// $Id: syssolver_mdagm.h,v 3.2 2009-06-02 15:56:40 bjoo Exp $
/*! \file
 *  \brief Disambiguator for MdagM system solvers
 */

#ifndef __syssolver_mdagm_h__
#define __syssolver_mdagm_h__

#include "linearop.h"
#include "handle.h"
#include "syssolver.h"
#include "update/molecdyn/predictor/chrono_predictor.h"

namespace Chroma
{
  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct MdagMSystemSolver : virtual public SystemSolver<T>
  {    

    virtual SystemSolverResults_t operator() (T& psi, const T& chi) const = 0;

    //! Return the subset on which the operator acts
    virtual const Subset& subset() const = 0;
    virtual SystemSolverResults_t operator()(T& psi, 
					     const T& chi,
					     AbsChronologicalPredictor4D<T>& predictor) const = 0;

  };


  //! SystemSolver disambiguator
  /*! This struct is solely to disambiguate the type of SystemSolvers */
  template<typename T>
  struct MdagMSystemSolverArray : virtual public SystemSolverArray<T>
  {
  };

}


#endif
