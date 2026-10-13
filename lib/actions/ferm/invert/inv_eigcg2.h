// -*- C++ -*-
/*! \file
 *  \brief Conjugate-Gradient algorithm with eigenstd::vector acceleration
 */

#ifndef __inv_eig_cg2_h__
#define __inv_eig_cg2_h__

#include "linearop.h"
#include "syssolver.h"
#include "actions/ferm/invert/containers.h"

#if ! defined (QDP_IS_QDPJIT2)

namespace Chroma 
{

  //! Conjugate-Gradient (CGNE) with eigenstd::vector acceleration
  /*! \ingroup invert
   * @{
   */
  namespace InvEigCG2Env
  {
    // LatticeFermionF
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperator<LatticeFermionF>& A,
			const multi1d<LatticeFermionF>& evec,
			int Nvecs);
    
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperator<LatticeFermionF>& A,
			const multi1d<LatticeFermionF>& evec,
			const multi1d<Double>& eval,
			int Nvecs,int NgoodEvecs) ;
    
    SystemSolverResults_t InvEigCG2(const LinearOperator<LatticeFermionF>& A,
				    LatticeFermionF& x, 
				    const LatticeFermionF& b,
				    multi1d<Double>& eval, 
				    multi1d<LatticeFermionF>& evec,
				    int Neig,
				    int Nmax,
				    const Real& RsdCG, int MaxCG, const int t);
  
    SystemSolverResults_t vecPrecondCG(const LinearOperator<LatticeFermionF>& A, 
				       LatticeFermionF& x, 
				       const LatticeFermionF& b, 
				       const multi1d<Double>& eval, 
				       const multi1d<LatticeFermionF>& evec, 
				       int startV, int endV,
				       const Real& RsdCG, int MaxCG);

    void InitGuess(const LinearOperator<LatticeFermionF>& A, 
		   LatticeFermionF& x, 
		   const LatticeFermionF& b, 
		   const multi1d<Double>& eval, 
		   const multi1d<LatticeFermionF>& evec, 
		   int& n_count);
  
    void InitGuess(const LinearOperator<LatticeFermionF>& A, 
		   LatticeFermionF& x, 
		   const LatticeFermionF& b, 
		   const multi1d<Double>& eval, 
		   const multi1d<LatticeFermionF>& evec, 
		   int N, // number of vectors to use
		   int& n_count);


    // LatticeFermionD
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperator<LatticeFermionD>& A,
			const multi1d<LatticeFermionD>& evec,
			int Nvecs);

    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperator<LatticeFermionD>& A,
			const multi1d<LatticeFermionD>& evec,
			const multi1d<Double>& eval,
			int Nvecs,int NgoodEvecs) ;

    SystemSolverResults_t InvEigCG2(const LinearOperator<LatticeFermionD>& A,
				    LatticeFermionD& x, 
				    const LatticeFermionD& b,
				    multi1d<Double>& eval, 
				    multi1d<LatticeFermionD>& evec,
				    int Neig,
				    int Nmax,
				    const Real& RsdCG, int MaxCG, const int f);
  
    SystemSolverResults_t vecPrecondCG(const LinearOperator<LatticeFermionD>& A, 
				       LatticeFermionD& x, 
				       const LatticeFermionD& b, 
				       const multi1d<Double>& eval, 
				       const multi1d<LatticeFermionD>& evec, 
				       int startV, int endV,
				       const Real& RsdCG, int MaxCG);

    void InitGuess(const LinearOperator<LatticeFermionD>& A, 
		   LatticeFermionD& x, 
		   const LatticeFermionD& b, 
		   const multi1d<Double>& eval, 
		   const multi1d<LatticeFermionD>& evec, 
		   int& n_count);
  
    void InitGuess(const LinearOperator<LatticeFermionD>& A, 
		   LatticeFermionD& x, 
		   const LatticeFermionD& b, 
		   const multi1d<Double>& eval, 
		   const multi1d<LatticeFermionD>& evec, 
		   int N, // number of vectors to use
		   int& n_count);

  } // namespace EigCG2Env

  /*! @} */  // end of group invert

  
}// End Namespace Chroma

#endif 
#endif 
