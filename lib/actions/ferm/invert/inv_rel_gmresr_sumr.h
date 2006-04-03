// -*- C++ -*-
// $Id: inv_rel_gmresr_sumr.h,v 3.0 2006-04-03 04:58:49 edwards Exp $
/*! \file
 *  \brief Relaxed GMRESR algorithm of the Wuppertal Group
 */

#ifndef __inv_rel_gmresr_sumr__
#define __inv_rel_gmresr_sumr__

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvRelGMRESR_SUMR(const LinearOperator<T>& PrecU,
		       const Complex& zeta,
		       const Real& rho,
		       const LinearOperator<T>& UnprecU,
		       const T& b,
		       T& x,
		       const Real& epsilon, 
		       const Real& epsilon_prec,
		       int MaxGMRESR, 
		       int MaxGMRESRPrec,
		       int& n_count);

}  // end namespace Chroma


#endif
