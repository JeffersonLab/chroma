// -*- C++ -*-
// $Id: inv_gmresr_cg_array.h,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Relaxed GMRESR algorithm of the Wuppertal Group
 */

#ifndef __inv_gmresr_cg_array_
#define __inv_gmresr_cg_array_

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvGMRESR_CG(const LinearOperatorArray<T>& PrecMM,
		  const LinearOperatorArray<T>& PrecM,
		  const LinearOperatorArray<T>& UnprecM,
		  const multi1d<T>& b,
		  multi1d<T>& x,
		  const Real& epsilon, 
		  const Real& epsilon_prec,
		  int MaxGMRESR, 
		  int MaxGMRESRPrec,
		  int& n_count);

}  // end namespace Chroma

#endif
