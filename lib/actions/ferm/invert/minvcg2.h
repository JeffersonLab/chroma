// -*- C++ -*-
// $Id: minvcg2.h,v 3.1 2008-03-10 17:32:40 bjoo Exp $
/*! \file
 *  \brief Multishift Conjugate-Gradient algorithm for a Linear Operator
 */

#ifndef MINVCG2_INCLUDE
#define MINVCG2_INCLUDE

#include "linearop.h"

namespace Chroma 
{


  /*! \ingroup invert */
  template<typename T>
  void MInvCG2(const LinearOperator<T>& M, 
	      const T& chi, 
	      multi1d<T>& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);


  void MInvCG2(const LinearOperator<LatticeFermionF>& M, 
	      const LatticeFermionF& chi, 
	      multi1d<LatticeFermionF>& psi,
	      const multi1d<RealF>& shifts, 
	      const multi1d<RealF>& RsdCG,
	      int MaxCG,
	      int &n_count);


  void MInvCG2(const LinearOperator<LatticeFermionD>& M, 
	      const LatticeFermionD& chi, 
	      multi1d<LatticeFermionD>& psi,
	      const multi1d<RealD>& shifts, 
	      const multi1d<RealD>& RsdCG,
	      int MaxCG,
	      int &n_count);

  /*! \ingroup invert */
  template<typename T, typename P, typename Q>
  void MInvCG2(const DiffLinearOperator<T,P,Q>& M,
	      const T& chi, 
	      multi1d<T>& psi,
	      const multi1d<Real>& shifts, 
	      const multi1d<Real>& RsdCG,
	      int MaxCG,
	      int &n_count);

  void MInvCG2(const DiffLinearOperator<LatticeFermionF,multi1d<LatticeColorMatrixF>,multi1d<LatticeColorMatrixF> >& M,
	      const LatticeFermionF& chi, 
	      multi1d<LatticeFermionF>& psi,
	      const multi1d<RealF>& shifts, 
	      const multi1d<RealF>& RsdCG,
	      int MaxCG,
	      int &n_count);

  void MInvCG2(const DiffLinearOperator<LatticeFermionD,multi1d<LatticeColorMatrixD>,multi1d<LatticeColorMatrixD> >& M,
	      const LatticeFermionD& chi, 
	      multi1d<LatticeFermionD>& psi,
	      const multi1d<RealD>& shifts, 
	      const multi1d<RealD>& RsdCG,
	      int MaxCG,
	      int &n_count);

}  // end namespace Chroma


#endif
