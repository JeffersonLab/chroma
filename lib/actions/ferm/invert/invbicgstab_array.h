// -*- C++ -*-
// $Id: invbicgstab_array.h,v 1.2 2005-01-14 20:13:05 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab_array__
#define __invbicgstab_array__

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvBiCGStab(const LinearOperator< multi1d<T> >& A,
		 const multi1d<T>& chi,
		 multi1d<T>& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count);

}  // end namespace Chroma

#endif
