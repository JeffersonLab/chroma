// -*- C++ -*-
// $Id: lapack_stub.cc,v 1.8 2009-08-03 15:19:19 bjoo Exp $
/*! \file
 *  \brief Stubs for QDP interface to Lapack lib
 */

#include "qdp-lapack.h"

using namespace QDP;

namespace QDPLapack
{
  int zheev(char& jobz,
	    char& uplo,
	    const int& N,
	    multi2d<DComplex>& A,
	    //const int& lda,
	    multi1d<Double>& w)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N,
	    multi2d<DComplex>& A,
	    //const int& lda,
	    multi1d<Double>& w,
	    multi1d<DComplex>& Work,
	    //const int& LWork,
	    multi1d<Double>& RWork)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,
	    multi1d<Double>& w)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zgeqrf(const int M,
	     const int N,
	     multi2d<DComplex>& A,
	     multi1d<DComplex>& TAU)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zunmqr(char& side,
	     char& trans,
	     const int M,
	     const int N,
	     multi2d<DComplex>& A,
	     multi1d<DComplex>& TAU,
	     multi2d<DComplex>& C)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
  int zunmqr2(char& side,
	      char& trans,
	      const int M,
	      const int N,
	      const int K,
	      multi2d<Complex64>& A, //input
	      multi1d<Complex64>& TAU, // some strange LAPACK beast
	      multi2d<Complex64>& C //input/output
	      )
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
  int zunmqrv(char& side,
	      char& trans,
	      const int M,
	      const int K,
	      multi2d<Complex64>& A, //input
	      multi1d<Complex64>& TAU, // some strange LAPACK beast
	      multi1d<Complex64>& C //input/output
	      )
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
  int zungqr(const int M, // The number of rows in Q
	     const int N,  // The number of columns in Q
	     const int K, // The number of reflections in Tau
	     multi2d<DComplex>& A, //input MxN matrix
	     multi1d<DComplex>& TAU // some strange LAPACK beast
	     )
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
  int zpotrf(char &uplo, int N,  multi2d<DComplex>& A, int LDA, int& info){
    
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;

  }
  int zpotrf(char &uplo, multi2d<DComplex>& A, int& info){
    int LDA = A.size1() ; // Assumes square matrix
    
    return zpotrf(uplo,LDA,A,LDA,info);
    
  }


  int cpotrf(char &uplo, int N,  multi2d<Complex>& A, int LDA, int& info){
    
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;

  }
  int cpotrf(char &uplo, multi2d<Complex>& A, int& info){
    int LDA = A.size1() ; // Assumes square matrix
    
    return cpotrf(uplo,LDA,A,LDA,info);
    
  }
  

  
  int zpotrs(char &uplo, int N, int nrhs,  
	     multi2d<DComplex>& A, int LDA, 
	     multi2d<DComplex>& B, int LDB, int& info){

    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
  int zpotrs(char &uplo, multi2d<DComplex>& A,  multi2d<DComplex>& B,  
	     int& info){
    
    int LDA = A.size1() ; // Assumes square matrix
    int LDB = B.size1() ; // second index is fast and its size is size1
    int nrhs = B.size2() ;
    
    return zpotrs(uplo,LDA,nrhs,A,LDA,B,LDB,info);
  }

 int zpotrs(char &uplo, int N, int nrhs,  
	     multi1d<Complex64>& A, int LDA, 
	     multi2d<Complex64>& B, int LDB, int& info){
   QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
 }

  int cpotrs(char &uplo, int N, int nrhs,  
	     multi1d<Complex>& A, int LDA, 
	     multi2d<Complex>& B, int LDB, int& info){

    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  
  }
  
  int cpotrs(char &uplo, multi1d<Complex>& A,  multi2d<Complex>& B,  
	     int& info){
    
    int LDA = A.size1() ; // Assumes square matrix
    int LDB = B.size1() ; // second index is fast and its size is size1
    int nrhs = B.size2() ;
    
    return cpotrs(uplo,LDA,nrhs,A,LDA,B,LDB,info);
    
  }
  
  
  void dsteqr(char *, int *, double *, double *, double *, int *,
               double *, int *)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__); 
  }

  void dlartg(Real64& F, Real64& G, Real64& CS, Real64& SN, Real64& R)
  {
    // Stub function
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
  }
  
    // ZGETRF
  int zgetrf(int M, int N, multi2d<Complex64>& A, int LDA,
	      multi1d<int>& ipiv, int& info)
  {
    // Stub function
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  void zgetrs(char &trans, int n, int nrhs,
	     multi2d<Complex64>& A, int lda, 
	     multi1d<int>& ipiv,  
	     multi1d<Complex64>& B,
	     int ldb, int& info) 
  {
    // Stub function
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
  }

  int zgeev(int n,
	    multi2d<DComplex>& A, 
	    multi1d<DComplex>& evals,
	    multi2d<DComplex>& evecs)
  {
    // Stub function
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  // The following are BLAS stubs:

  int zgemm(char &transa,
	    char &transb,
	    int m,
	    int n,
	    int k,
	    Complex64 alpha,
	    multi2d<Complex64>& A,  // input
	    int lda,
	    multi2d<Complex64>& B,  // input
	    int ldb,
	    Complex64 beta,
	    multi2d<Complex64>& C,  //input/output
	    int ldc   
    )
  {
    // Stub function
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zgemv(char &trans, 
	    const int& m, const int& n,
	    multi2d<Complex64>& A, 
	    int lda,
	    multi1d<Complex64>& x,
	    multi1d<Complex64>& y
	    )
  {
    // Stub function
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }


} // namespace Lapack


