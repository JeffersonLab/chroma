// -*- C++ -*-
// $Id: invbicgstab.h,v 1.2 2005-01-14 20:13:05 edwards Exp $
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
