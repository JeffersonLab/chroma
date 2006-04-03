#ifndef INV_MINRES_ARRAY_H
#define INV_MINRES_ARRAY_H

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvMINRES(const LinearOperatorArray<T>& A,
	       const multi1d<T>& chi,
	       multi1d<T>& psi,
	       const Real& RsdCG,
	       int MaxCG,
	       int& n_count);

}  // end namespace Chroma

#endif
