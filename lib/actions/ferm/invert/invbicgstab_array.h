// -*- C++ -*-
// $Id: invbicgstab_array.h,v 3.0 2006-04-03 04:58:49 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab_array__
#define __invbicgstab_array__

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvBiCGStab(const LinearOperatorArray<T>& A,
		 const multi1d<T>& chi,
		 multi1d<T>& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count);

}  // end namespace Chroma

#endif
