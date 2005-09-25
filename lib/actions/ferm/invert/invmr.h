// -*- C++ -*-
// $Id: invmr.h,v 2.0 2005-09-25 21:04:27 edwards Exp $

#ifndef __invmr_h__
#define __invmr_h__

#include "linearop.h"

namespace Chroma {

void InvMR(const LinearOperator& A,
	   const LatticeFermion& chi,
	   LatticeFermion& psi,
	   const Real& RsdMR, 
	   const Real& MRovpar;
	   int MaxCG, 
	   int& n_count);

}  // end namespace Chroma

#endif
