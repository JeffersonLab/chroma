// -*- C++ -*-
// $Id: invbicgstab_array.h,v 2.0 2005-09-25 21:04:27 edwards Exp $
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
