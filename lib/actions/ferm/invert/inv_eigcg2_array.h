// -*- C++ -*-
// $Id: inv_eigcg2_array.h,v 1.1 2008-01-13 21:08:13 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm with eigenvector acceleration
 */

#ifndef __inv_eigcg2_array_h__
#define __inv_eigcg2_array_h__

#include "linearop.h"
#include "syssolver.h"
#include "actions/ferm/invert/containers.h"

namespace Chroma 
{

  //! Conjugate-Gradient (CGNE) with eigenvector acceleration
  /*! \ingroup invert
   * @{
   */
  namespace InvEigCG2ArrayEnv
  {
    // LatticeFermionF
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperatorArray<LatticeFermionF>& A,
			const multi2d<LatticeFermionF>& evec,
			int Nvecs);

    SystemSolverResults_t InvEigCG2(const LinearOperatorArray<LatticeFermionF>& A,
				    multi1d<LatticeFermionF>& x, 
				    const multi1d<LatticeFermionF>& b,
				    multi1d<Double>& eval, 
				    multi2d<LatticeFermionF>& evec,
				    int Neig,
				    int Nmax,
				    const Real& RsdCG, int MaxCG);
  
    SystemSolverResults_t vecPrecondCG(const LinearOperatorArray<LatticeFermionF>& A, 
				       multi1d<LatticeFermionF>& x, 
				       const multi1d<LatticeFermionF>& b, 
				       const multi1d<Double>& eval, 
				       const multi2d<LatticeFermionF>& evec, 
				       int startV, int endV,
				       const Real& RsdCG, int MaxCG);

    void InitGuess(const LinearOperatorArray<LatticeFermionF>& A, 
		   multi1d<LatticeFermionF>& x, 
		   const multi1d<LatticeFermionF>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionF>& evec, 
		   int& n_count);
  
    void InitGuess(const LinearOperatorArray<LatticeFermionF>& A, 
		   multi1d<LatticeFermionF>& x, 
		   const multi1d<LatticeFermionF>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionF>& evec, 
		   int N, // number of vectors to use
		   int& n_count);


    // multi1d<LatticeFermionD>
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperatorArray<LatticeFermionD>& A,
			const multi2d<LatticeFermionD>& evec,
			int Nvecs);

    SystemSolverResults_t InvEigCG2(const LinearOperatorArray<LatticeFermionD>& A,
				    multi1d<LatticeFermionD>& x, 
				    const multi1d<LatticeFermionD>& b,
				    multi1d<Double>& eval, 
				    multi2d<LatticeFermionD>& evec,
				    int Neig,
				    int Nmax,
				    const Real& RsdCG, int MaxCG);
  
    SystemSolverResults_t vecPrecondCG(const LinearOperatorArray<LatticeFermionD>& A, 
				       multi1d<LatticeFermionD>& x, 
				       const multi1d<LatticeFermionD>& b, 
				       const multi1d<Double>& eval, 
				       const multi2d<LatticeFermionD>& evec, 
				       int startV, int endV,
				       const Real& RsdCG, int MaxCG);

    void InitGuess(const LinearOperatorArray<LatticeFermionD>& A, 
		   multi1d<LatticeFermionD>& x, 
		   const multi1d<LatticeFermionD>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionD>& evec, 
		   int& n_count);
  
    void InitGuess(const LinearOperatorArray<LatticeFermionD>& A, 
		   multi1d<LatticeFermionD>& x, 
		   const multi1d<LatticeFermionD>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionD>& evec, 
		   int N, // number of vectors to use
		   int& n_count);

  } // namespace EigCG2Env

  /*! @} */  // end of group invert

  
}// End Namespace Chroma

#endif 
