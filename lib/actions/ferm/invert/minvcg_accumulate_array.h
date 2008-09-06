// -*- C++ -*-
// $Id: minvcg_accumulate_array.h,v 3.1 2008-09-06 18:35:35 bjoo Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#ifndef __minvcg_accumulate_array_h__
#define __minvcg_accumulate_array_h__

#include "linearop.h"


using namespace QDP;

namespace Chroma 
{

  /*! \ingroup invert */
  template<typename T>
  void MInvCGAccum(const LinearOperatorArray<T>& A,
		   const multi1d<T>& chi, 
		   multi1d<T>& psi,
		   const Real& norm,
		   const multi1d<Real>& residues,
		   const multi1d<Real>& poles, 
		   const Real& RsdCG,
		   const int MaxCG,
		   int &n_count);


  /*! \ingroup invert */
  template<typename T, typename P, typename Q>
  void MInvCGAccum(const DiffLinearOperatorArray<T,P,Q>& A,
		   const multi1d<T>& chi, 
		   multi1d<T>& psi,
		   const Real& norm,
		   const multi1d<Real>& residues,
		   const multi1d<Real>& poles, 
		   const Real& RsdCG,
		   const int MaxCG,
		   int &n_template);

  
}  // end namespace Chroma


#endif
