// -*- C++ -*-
// $Id: lapack_wrapper.cc,v 1.1 2007-09-26 02:50:35 kostas Exp $

#include "lapack_wrapper.h"


// NEEDS A LOT OF CLEAN UP
namespace Lapack
{
 
  // hide lapack's temporaries...
  int zheev(char& jobz,
	    char& uplo,
	    const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w){
    /**
    char c_jobz = *jobz;
    char c_uplo = *uplo;
    **/

    int LWork = 2*N-1;
    multi1d<DComplex> Work(LWork);
    multi1d<Double> RWork(3*N-2);
    
   int lda = A.size1() ; 
        

    w.resize(N) ;

    int info;
    int r = zheev_(&jobz,&uplo,
		   (int *)&N,
		   &A(0,0),
		   &lda,
		   &w[0],
		   &Work[0],
		   &LWork,
		   &RWork[0],
		   &info) ;
    
    if(info){
      QDPIO::cerr<<"Lapack::zheev returned with exit code: "<<info<<endl ;
      exit(1) ;
    }

    return r ;
}

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w,
	    multi1d<DComplex>& Work, // Should be length LWork >= max(1,2*N-1)
	    //const int& LWork,
	    multi1d<Double>& RWork // Should be length max(1,3*N-2)
	    ){
    /**
    char c_jobz = *jobz;
    char c_uplo = *uplo;
    **/
    int N = A.size1();
    int lda = N ; 
    int LWork = Work.size();
        
    w.resize(N) ;

    int info ;
    int r =  zheev_(&jobz,&uplo,
		    &N,
		    &A(0,0),
		    &lda,
		    &w[0],
		    &Work[0],
		    &LWork,
		    &RWork[0],
		    &info) ;
    
    if(info){
      QDPIO::cerr<<"Lapack::zheev returned with exit code: "<<info<<endl ;
      exit(1) ;
    }

    return r ;

  }

  // hide lapack's temporaries...
  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w){

    int N = A.size1();
    int LWork = 2*N-1;
    multi1d<DComplex> Work(LWork);
    multi1d<Double> RWork(3*N-2);
    
    return zheev(jobz,uplo,A,w,Work,RWork) ;
  }

  int zgeqrf(const int M, // The vector length
	     const int N, // The number of vectors
	     multi2d<DComplex>& A, // the array containing the vectors
	     multi1d<DComplex>& TAU // some strange LAPACK beast
	     )
  {
    if(N>M)
      TAU.resize(M);
    else
      TAU.resize(N);
	
    int LWork = N ;
    multi1d<DComplex> Work(LWork);
    

    int lda = A.size1(); // need to check which is LDA size1 or size2
    //But sice I am using square matrices this is OK
    int info ;
    int r = zgeqrf_((int *)&M,
		    (int *)&N,
		    &A(0,0),
		    &lda,
		    &TAU[0],
		    &Work[0],
		    &LWork,
		    &info) ;

    
    if(info){
      QDPIO::cerr<<"Lapack::zgeqrf returned with exit code: "<<info<<endl ;
      exit(1) ;
    }

    return r ;

  }

 /** 
   *  Purpose
   *  =======
   *
   *  ZUNMQR overwrites the general complex M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'C':      Q**H * C       C * Q**H
   *
   *  where Q is a complex unitary matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   **/
  int zunmqr(char& side,
	     char& trans,
	     const int M,
	     const int N,
	     multi2d<DComplex>& A, //input
	     multi1d<DComplex>& TAU, // some strange LAPACK beast
	     multi2d<DComplex>& C //input/output
	     ){
    /**
    char c_side = *side ;
    char c_trans = *trans ;
    **/

    int K = TAU.size();
    int lda = A.size1();
    int ldc = C.size1();
    int LWork = -1 ;
    if(side == 'R')
      LWork = M ;
    if(side == 'L')
      LWork = N ;

    /**
    cout<<"ZUNMQR->side  : "<<side<<endl  ;
    cout<<"ZUNMQR->trans : "<<trans<<endl ;
    cout<<"ZUNMQR->M     : "<<M<<endl ;
    cout<<"ZUNMQR->N     : "<<N<<endl ;
    cout<<"ZUNMQR->K     : "<<K<<endl ;
    cout<<"ZUNMQR->lda   : "<<lda<<endl ;
    cout<<"ZUNMQR->ldc   : "<<ldc<<endl ;
    cout<<"ZUNMQR->LWork : "<<LWork<<endl ;
    **/

    multi1d<DComplex> Work(LWork);
    
    int info ;
    int r = zunmqr_(&side, &trans,
		    (int *)&M, 
		    (int *)&N, 
		    (int *)&K, 
		    &A(0,0),
		    &lda, 
		    &TAU[0], 
		    &C(0,0),
		    &ldc,
		    &Work[0],
		    &LWork,
		    &info);

    if(info){
      QDPIO::cerr<<"Lapack::zunmqr returned with exit code: "<<info<<endl ;
      exit(1) ;
    }

    return r ;
  }
    
  
}// End Namespace Lapack


