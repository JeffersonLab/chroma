// -*- C++ -*-
// $Id: minvcg.h,v 1.5 2005-01-14 20:13:05 edwards Exp $

#ifndef MINVCG_INCLUDE
#define MINVCG_INCLUDE

#include "linearop.h"

namespace Chroma {


template<typename T>
void MInvCG(const LinearOperator<T>& A, 
	    const T& chi, 
	    multi1d<T>& psi,
	    const multi1d<Real>& shifts, 
	    const multi1d<Real>& RsdCG,
	    int MaxCG,
	    int &n_count);

}  // end namespace Chroma


#endif
