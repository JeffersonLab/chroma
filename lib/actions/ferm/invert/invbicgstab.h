// -*- C++ -*-
// $Id: invbicgstab.h,v 3.3 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab__
#define __invbicgstab__

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
  InvBiCGStab(const LinearOperator<T>& A,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      int MaxBiCGStab,
	      enum PlusMinus isign);

 

  /*! @} */  // end of group invert
	    
}  // end namespace Chroma

#endif
