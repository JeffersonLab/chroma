// -*- C++ -*-
// $Id: invcg2.h,v 1.1 2003-04-08 21:36:13 edwards Exp $

#ifndef __invcg2__
#define __invcg2__

#include "linearop.h"

void InvCG2(const LinearOperator& M,
	    const LatticeFermion& chi,
	    LatticeFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count);

#endif
