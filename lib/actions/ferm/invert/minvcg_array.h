// -*- C++ -*-
// $Id: minvcg_array.h,v 1.1 2005-01-28 02:14:28 edwards Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#ifndef __minvcg_array_h__
#define __minvcg_array_h__

#include "linearop.h"

namespace Chroma 
{

  template<typename T>
  void MInvCG(const LinearOperator< multi1d<T> >& A, 
	      const multi1d<T>& chi, 
	      multi1d< multi1d<T> >& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);

}  // end namespace Chroma


#endif
