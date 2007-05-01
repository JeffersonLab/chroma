// -*- C++ -*-
// $Id: invbicgstab_array.h,v 3.1 2007-05-01 12:47:24 bjoo Exp $
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
	    const Real& Relax,
	    int MaxBiCGStab);

}  // end namespace Chroma

#endif
