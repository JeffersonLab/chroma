// -*- C++ -*-
// $Id: invcg2_array.h,v 1.1 2003-11-13 04:13:53 edwards Exp $

#ifndef __invcg2_array__
#define __invcg2_array__

#include "linearop.h"

template<typename T>
void InvCG2(const LinearOperator< multi1d<T> >& M,
	    const multi1d<T>& chi,
	    multi1d<T>& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count);

#endif
