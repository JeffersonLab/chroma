// -*- C++ -*-
// $Id: inv_eigcg2.h,v 1.3 2007-10-09 18:07:14 edwards Exp $

#ifndef _INV_EIG_CG2_H
#define _INV_EIG_CG2_H

#include "linearop.h"
#include "syssolver.h"
#include "containers.h"


// NEEDS A LOT OF CLEAN UP
namespace Chroma 
{

  template<typename T>
  void SimpleGramSchmidt(multi1d<T>& vec, 
			 int f,
			 int t,
			 const Subset& sub)
  {
    for(int i(0);i<f;i++){// normalize the first vectors...
      vec[i][sub] /= Real(sqrt(norm2(vec[i],sub))) ;
    }
    if(!(t<=vec.size())){
      QDPIO::cerr<<"SimpleGramSchmidt:: f="<<f<<" t="<<t<<" vec.size()="<<vec.size()<<endl;
      QDPIO::cerr<<"SimpleGramSchmidt:: Out of bound!\n";
      exit(1);
    }
    for(int i(f);i<t;i++)
    { // now orthonormalize ther rest
      for(int k(0);k<i;k++){
	DComplex dcc = innerProduct(vec[k], vec[i], sub);
	Complex cc = dcc ;
	//cout<<"GramS: "<<cc<<" "<<dcc<<endl;
	vec[i][sub] = vec[i]  - cc*vec[k] ;
      }
      Double in = 1.0/sqrt(norm2(vec[i],sub)) ;
      vec[i][sub] *= Real(in); 
    }
  }

  void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
		      const LinearOperator<LatticeFermion>& A,
		      const multi1d<LatticeFermion>& evec,
		      int Nvecs) ;

  SystemSolverResults_t InvEigCG2(const LinearOperator<LatticeFermion>& A,
				  LatticeFermion& x, 
				  const LatticeFermion& b,
				  multi1d<Double>& eval, 
				  multi1d<LatticeFermion>& evec,
				  int Neig,
				  int Nmax,
				  const Real& RsdCG, int MaxCG) ;
  
  SystemSolverResults_t vecPrecondCG(const LinearOperator<LatticeFermion>& A, 
				     LatticeFermion& x, 
				     const LatticeFermion& b, 
				     const multi1d<Double>& eval, 
				     const multi1d<LatticeFermion>& evec, 
				     int startV, int endV,
				     const Real& RsdCG, int MaxCG) ;

  void InitGuess(const LinearOperator<LatticeFermion>& A, 
		 LatticeFermion& x, 
		 const LatticeFermion& b, 
		 const multi1d<Double>& eval, 
		 const multi1d<LatticeFermion>& evec, 
		 int& n_count) ;
  
  void InitGuess(const LinearOperator<LatticeFermion>& A, 
		 LatticeFermion& x, 
		 const LatticeFermion& b, 
		 const multi1d<Double>& eval, 
		 const multi1d<LatticeFermion>& evec, 
		 int N, // number of vectors to use
		 int& n_count) ;

 void InitCG(const LinearOperator<LatticeFermion>& A, 
	     LatticeFermion& x, 
	     const LatticeFermion& b, 
	     const multi1d<Double>& eval, 
	     const multi1d<LatticeFermion>& evec, 
	     const Real& RsdCG, int MaxCG, int& n_count) ;
  

  
}// End Namespace Chroma

#endif 
