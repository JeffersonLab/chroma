// -*- C++ -*-
// $Id: inv_rel_gmresr_cg.h,v 1.2 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief Relaxed GMRESR algorithm of the Wuppertal Group
 */

#ifndef __inv_rel_gmresr_cg__
#define __inv_rel_gmresr_cg__

#include "chromabase.h"
#include "linearop.h"


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

#endif
