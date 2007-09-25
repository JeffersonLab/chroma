// -*- C++ -*-
// $Id: inv_eigcg2.h,v 1.1 2007-09-25 21:17:12 edwards Exp $

#ifndef _INV_EIG_CG2_H
#define _INV_EIG_CG2_H

#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
#include "lmdagm.h"
#include "containers.h"
#include "simpleGramSchmidt.h"
#include "lapack_wrapper.h"


// NEEDS A LOT OF CLEAN UP
namespace Chroma 
{

  void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
		      LinearOperator<LatticeFermion>& A,
		      const multi1d<LatticeFermion>& evec,
		      const int Nvecs) ;


  void InvEigCG2(LinearOperator<LatticeFermion>& A,
		 LatticeFermion& x, 
		 const LatticeFermion& b,
		 multi1d<Double>& eval, 
		 multi1d<LatticeFermion>& evec,
		 const int Neig,
		 const int Nmax,
		 const Real RsdCG, const int MaxCG,int& n_count) ;

  void vecPrecondCG(LinearOperator<LatticeFermion>& A, 
		    LatticeFermion& x, 
		    const LatticeFermion& b, 
		    const multi1d<Double>& eval, 
		    const multi1d<LatticeFermion>& evec, 
		    const Real RsdCG, const int MaxCG, int& n_count) ;

  void InitGuess(LinearOperator<LatticeFermion>& A, 
		 LatticeFermion& x, 
		 const LatticeFermion& b, 
		 const multi1d<Double>& eval, 
		 const multi1d<LatticeFermion>& evec, 
		 int& n_count) ;
  
  void InitGuess(LinearOperator<LatticeFermion>& A, 
		 LatticeFermion& x, 
		 const LatticeFermion& b, 
		 const multi1d<Double>& eval, 
		 const multi1d<LatticeFermion>& evec, 
		 const int N, // number of vectors to use
		 int& n_count) ;

 void InitCG(LinearOperator<LatticeFermion>& A, 
	     LatticeFermion& x, 
	     const LatticeFermion& b, 
	     const multi1d<Double>& eval, 
	     const multi1d<LatticeFermion>& evec, 
	     const Real RsdCG, const int MaxCG, int& n_count) ;
  

  
}// End Namespace Chroma

#endif 
