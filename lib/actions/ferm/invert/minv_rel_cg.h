// -*- C++ -*-
// $Id: minv_rel_cg.h,v 1.1 2004-05-21 15:31:50 bjoo Exp $

#ifndef MINV_REL_CG_INCLUDE
#define MINV_REL_CG_INCLUDE

#include "chromabase.h"
#include "linearop.h"

template<typename T>
void MInvRelCG(const ApproxLinearOperator<T>& A, 
	       const T& chi, 
	       multi1d<T>& psi,
	       const multi1d<Real>& shifts, 
	       const multi1d<Real>& RsdCG,
	       int MaxCG,
	       int &n_count);

#endif
