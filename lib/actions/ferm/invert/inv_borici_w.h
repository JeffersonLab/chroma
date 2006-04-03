#ifndef INV_BORICI_H
#define INV_BORICI_H

#include "linearop.h"

namespace Chroma {

template<typename T> 
void InvBorici( const LinearOperator<T>& D_4,
		const LinearOperatorArray<T>& D_5,
		const LinearOperatorArray<T>& D_dag_D_5,
		const T& b,
		T& x,
		const Real& tol,
		const Real& tol_1,
		const int MaxIter,
		const int MaxIter5D,
		const Real& m,
		int& n_iters);

}  // end namespace Chroma

#endif
