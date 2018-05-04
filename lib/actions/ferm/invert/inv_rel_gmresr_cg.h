// -*- C++ -*-
/*! \file
 *  \brief Relaxed GMRESR algorithm of the Wuppertal Group
 */

#ifndef __inv_rel_gmresr_cg__
#define __inv_rel_gmresr_cg__

#include "linearop.h"

namespace Chroma {

template<typename T>
void InvRelGMRESR_CG(const LinearOperator<T>& PrecMM,
		     const LinearOperator<T>& UnprecMM,
		     const T& b,
		     T& x,
		     const Real& epsilon, 
		     const Real& epsilon_prec,
		     int MaxGMRESR, 
		     int MaxGMRESRPrec,
		     int& n_count);

}  // end namespace Chroma

#endif
