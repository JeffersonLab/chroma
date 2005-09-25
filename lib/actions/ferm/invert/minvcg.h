// -*- C++ -*-
// $Id: minvcg.h,v 2.0 2005-09-25 21:04:28 edwards Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

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
