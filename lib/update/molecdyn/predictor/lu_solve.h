// -*- C++ -*-
// $Id: lu_solve.h,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 * \brief LU solver
 *
 * Predictors for HMC
 */

#ifndef LU_SOLVE_H
#define LU_SOLVE_H

#include "chromabase.h"

using namespace QDP;

namespace Chroma 
{

  //! Solve M a = b by LU decomposition with partial pivoting
  /*! @ingroup predictor */
  void LUSolve( multi1d<DComplex>& a, 
		const multi2d<DComplex>& M, 
		const multi1d<DComplex>& b );

};

#endif
