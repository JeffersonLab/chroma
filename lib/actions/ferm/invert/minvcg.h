// -*- C++ -*-
// $Id: minvcg.h,v 1.4 2004-01-08 17:11:01 bjoo Exp $

#ifndef MINVCG_INCLUDE
#define MINVCG_INCLUDE

#include "chromabase.h"
#include "linearop.h"

template<typename T>
void MInvCG(const LinearOperator<T>& A, 
	    const T& chi, 
	    multi1d<T>& psi,
	    const multi1d<Real>& shifts, 
	    const multi1d<Real>& RsdCG,
	    int MaxCG,
	    int &n_count);

#endif
