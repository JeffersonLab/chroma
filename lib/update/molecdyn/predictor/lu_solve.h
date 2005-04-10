// -*- C++ -*-
// $Id: lu_solve.h,v 1.2 2005-04-10 21:46:43 edwards Exp $
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
