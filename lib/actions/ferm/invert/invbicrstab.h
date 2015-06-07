// -*- C++ -*-
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicrstab__
#define __invbicrstab__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma 
{

  //! Bi-CG stabilized
  /*! \ingroup invert
   *
   * @{
   */
  template<typename T>
  SystemSolverResults_t
  InvBiCRStab(const LinearOperator<T>& A,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      int MaxBiCGStab,
	      enum PlusMinus isign);

 

  /*! @} */  // end of group invert
	    
}  // end namespace Chroma

#endif
