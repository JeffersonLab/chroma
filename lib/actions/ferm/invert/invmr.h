// -*- C++ -*-
// $Id: invmr.h,v 1.3 2005-01-14 20:13:05 edwards Exp $

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
