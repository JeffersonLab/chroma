// -*- C++ -*-
// $Id: minvcg_array.h,v 3.0 2006-04-03 04:58:49 edwards Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#ifndef __minvcg_array_h__
#define __minvcg_array_h__

#include "linearop.h"

namespace Chroma 
{

  /*! \ingroup invert */
  template<typename T>
  void MInvCG(const LinearOperatorArray<T>& A,
	      const multi1d<T>& chi, 
	      multi1d< multi1d<T> >& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);


  /*! \ingroup invert */
  template<typename T, typename P, typename Q>
  void MInvCG(const DiffLinearOperatorArray<T,P,Q>& A,
	      const multi1d<T>& chi, 
	      multi1d< multi1d<T> >& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);

}  // end namespace Chroma


#endif
