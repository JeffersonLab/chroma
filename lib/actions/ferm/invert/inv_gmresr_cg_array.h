// -*- C++ -*-
// $Id: inv_gmresr_cg_array.h,v 1.2 2005-01-14 20:13:04 edwards Exp $
/*! \file
 *  \brief Relaxed GMRESR algorithm of the Wuppertal Group
 */

#ifndef __inv_gmresr_cg_array_
#define __inv_gmresr_cg_array_

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvGMRESR_CG(const LinearOperator< multi1d<T> >& PrecMM,
		  const LinearOperator< multi1d<T> >& PrecM,
		  const LinearOperator< multi1d<T> >& UnprecM,
		  const multi1d<T>& b,
		  multi1d<T>& x,
		  const Real& epsilon, 
		  const Real& epsilon_prec,
		  int MaxGMRESR, 
		  int MaxGMRESRPrec,
		  int& n_count);

}  // end namespace Chroma

#endif
