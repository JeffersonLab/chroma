// -*- C++ -*-
// $Id: minvcg.h,v 1.6 2005-01-28 02:14:05 edwards Exp $
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
