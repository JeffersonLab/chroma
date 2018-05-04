// -*- C++ -*-

#ifndef MINV_REL_CG_INCLUDE
#define MINV_REL_CG_INCLUDE

#include "linearop.h"

namespace Chroma {

template<typename T>
void MInvRelCG(const LinearOperator<T>& A, 
	       const T& chi, 
	       multi1d<T>& psi,
	       const multi1d<Real>& shifts, 
	       const multi1d<Real>& RsdCG,
	       int MaxCG,
	       int &n_count);

}  // end namespace Chroma

#endif
