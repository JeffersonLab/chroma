// -*- C++ -*-
// $Id: minvcg2_accum.h,v 3.1 2008-09-02 20:10:18 bjoo Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#ifndef MINVCG2_ACCUM_INCLUDE
#define MINVCG2_ACCUM_INCLUDE

#include "linearop.h"

namespace Chroma 
{


  /*! \ingroup invert */
  template<typename T>
  void MInvCG2Accum(const LinearOperator<T>& M, 
		    const T& chi, 
		    T& psi,
		    const Real& norm, 
		    const multi1d<Real>& residues,
		    const multi1d<Real>& poles, 
		    const Real&  RsdCG,
		    int MaxCG,
		    int &n_count);

  /*! \ingroup invert */
  template<typename T, typename P, typename Q>
  void MInvCG2Accum(const DiffLinearOperator<T,P,Q>& M,
		    const T& chi, 
		    T& psi,
		    const Real& norm,
		    const multi1d<Real>& residues,
		    const multi1d<Real>& poles, 
		    const Real& RsdCG,
		    int MaxCG,
		    int &n_count);

}  // end namespace Chroma


#endif
