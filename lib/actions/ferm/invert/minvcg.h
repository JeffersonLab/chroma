// -*- C++ -*-
// $Id: minvcg.h,v 3.0 2006-04-03 04:58:49 edwards Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#ifndef MINVCG_INCLUDE
#define MINVCG_INCLUDE

#include "linearop.h"

namespace Chroma 
{


  /*! \ingroup invert */
  template<typename T>
  void MInvCG(const LinearOperator<T>& A, 
	      const T& chi, 
	      multi1d<T>& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);

  /*! \ingroup invert */
  template<typename T, typename P, typename Q>
  void MInvCG(const DiffLinearOperator<T,P,Q>& A,
	      const T& chi, 
	      multi1d<T>& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);

}  // end namespace Chroma


#endif
