// -*- C++ -*-
// $Id: invbicgstab.h,v 3.1 2007-05-01 12:47:24 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab__
#define __invbicgstab__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma {

template<typename T>
SystemSolverResults_t
InvBiCGStab(const LinearOperator<T>& A,
	    const T& chi,
	    T& psi,
	    const Real& RsdBiCGStab,
	    int MaxBiCGStab,
	    enum PlusMinus isign);

	    
}  // end namespace Chroma

#endif
