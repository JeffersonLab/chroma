// -*- C++ -*-
// $Id: invcg2.h,v 1.3 2003-10-20 20:31:50 edwards Exp $

#ifndef __invcg2__
#define __invcg2__

#include "linearop.h"

template<typename T>
void InvCG2(const LinearOperator<T>& M,
	    const T& chi,
	    T& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count);

#endif
