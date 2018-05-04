// -*- C++ -*-
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab_array__
#define __invbicgstab_array__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma {

template<typename T>
SystemSolverResults_t 
InvBiCGStab(const LinearOperatorArray<T>& A,
	    const multi1d<T>& chi,
	    multi1d<T>& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab);

}  // end namespace Chroma

#endif
