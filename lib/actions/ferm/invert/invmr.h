// -*- C++ -*-
// $Id: invmr.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

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
