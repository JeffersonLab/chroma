// -*- C++ -*-
// $Id: inv_gmresr_cg_array.h,v 1.1 2004-05-27 11:37:44 bjoo Exp $
/*! \file
 *  \brief Relaxed GMRESR algorithm of the Wuppertal Group
 */

#ifndef __inv_gmresr_cg_array_
#define __inv_gmresr_cg_array_

#include "chromabase.h"
#include "linearop.h"


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

#endif
