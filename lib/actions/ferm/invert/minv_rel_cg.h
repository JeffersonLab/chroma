// -*- C++ -*-
// $Id: minv_rel_cg.h,v 1.2 2004-12-12 21:22:15 edwards Exp $

#ifndef MINV_REL_CG_INCLUDE
#define MINV_REL_CG_INCLUDE

#include "chromabase.h"
#include "linearop.h"

template<typename T>
void MInvRelCG(const LinearOperator<T>& A, 
	       const T& chi, 
	       multi1d<T>& psi,
	       const multi1d<Real>& shifts, 
	       const multi1d<Real>& RsdCG,
	       int MaxCG,
	       int &n_count);

#endif
