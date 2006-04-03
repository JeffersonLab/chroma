// -*- C++ -*-
// $Id: invbicgstab.h,v 3.0 2006-04-03 04:58:49 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab__
#define __invbicgstab__

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvBiCGStab(const LinearOperator<T>& A,
		 const T& chi,
		 T& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count);

}  // end namespace Chroma

#endif
