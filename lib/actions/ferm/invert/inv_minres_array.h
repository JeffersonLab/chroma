#ifndef INV_MINRES_ARRAY_H
#define INV_MINRES_ARRAY_H

#include "linearop.h"

template<typename T>
void InvMINRES(const LinearOperator< multi1d<T> >& A,
	       const multi1d<T>& chi,
	       multi1d<T>& psi,
	       const Real& RsdCG,
	       int MaxCG,
	       int& n_count);

#endif
