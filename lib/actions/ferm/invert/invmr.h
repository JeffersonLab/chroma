// -*- C++ -*-
// $Id: invmr.h,v 1.1 2003-04-14 16:35:01 edwards Exp $

#ifndef __invmr_h__
#define __invmr_h__

#include "linearop.h"

void InvMR(const LinearOperator& A,
	   const LatticeFermion& chi,
	   LatticeFermion& psi,
	   const Real& RsdMR, 
	   const Real& MRovpar;
	   int MaxCG, 
	   int& n_count);

#endif
