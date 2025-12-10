// -*- C++ -*-
// $Id: fortran_lapack.h,v 1.8 2009-07-31 14:09:31 bjoo Exp $
/*! \file
 *  \brief QDP interface to Lapack lib
 */
//MAC OS X
//#include <vecLib/vBLAS.h>
//#include <veclib/clapack.h>

// May need other includes for other OS 
#include "qdp.h"
using namespace QDP;

extern "C" {

  //LAPACK routines
  int zheev_(char *,char *,
	     int *, Complex64 *, int *, Real64 *, Complex64 *, 
	     int *, Real64 *, int *);
  
  int  zgeqrf_(int *, int *,
	       Complex64 *, int *,
	       Complex64 *, Complex64 *,
	       int *, int *) ;
  
  int zunmqr_(char *, char *,
	      int *, int *,
	      int *, Complex64 *,
	      int *, Complex64 *,  Complex64 *,
	      int *, Complex64 *,
	      int *, int *);

  int zungqr_(const int*,        // M
	      const int*,        // N
	      const int*,        // K
	      Complex64*,  // A
	      const int*,        // LDA
	      Complex64*,  // TAU
	      Complex64*,  // WORK
	      const int*,        // LWORK
	      int*);       // INFO

  int zhetrf_(char *uplo, int *n, Complex64 *a, 
	      int *lda, int *ipiv, Complex64 *work, int *lwork, 
	      int *info);

  int zgetrf_(int* M, int *N, Complex64* a, int* lda, int *ipiv, int *info);

  int zhetrs_(char *uplo, int *n, int *nrhs, 
	      Complex64 *a, int *lda, int *ipiv, Complex64 *b, 
	      int *ldb, int *info);

  void zgetrs_(char *TRANS,
	      int* N,
	      int* nrhs,
	      Complex64* A,
	      int* lda,
	      int* ipiv,
	      Complex64* B,
	      int* ldb,
	      int* info);

  // Cholesky factorization
  int zpotrf_(char *uplo, int *n, Complex64 *a, int *lda, int *info);
  int cpotrf_(char *uplo, int *n, Complex32 *a, int *lda, int *info);

  // Cholesky inverse
  int zpotri_(char *uplo, int *n, Complex64 *a, int *lda, int *info);
  int cpotri_(char *uplo, int *n, Complex32 *a, int *lda, int *info);

  // Cholesky linear system solve
  int zpotrs_(char *uplo, int *n, int *nrhs, Complex64 *a, int *lda, 
	      Complex64 *b, int *ldb, int *info);
  int cpotrs_(char *uplo, int *n, int *nrhs, Complex32 *a, int *lda, 
	      Complex32 *b, int *ldb, int *info);

  

  //BLAS routines  
  int cgemm_(char *TRANSA, char *TRANSB,int *M,int *N,int *K, 
             Complex32 *ALPHA,
             Complex32 *A, int *LDA, 
             Complex32 *B, int *LDB,
             Complex32 *BETA,
             Complex32 *C, int *LDC) ;

  int zgemm_(char *TRANSA, char *TRANSB,int *M,int *N,int *K, 
             Complex64 *ALPHA,
             Complex64 *A, int *LDA, 
             Complex64 *B, int *LDB,
             Complex64 *BETA,
             Complex64 *C, int *LDC) ;

  int cgemv_(char *trans, int *M, int *N, Complex32 *alpha, 
	       Complex32 *A, int *lda, Complex32 *x, int *incx, 
	       Complex32 *beta, Complex32 *y, int *incy);

  int zgemv_(char *trans, int *M, int *N, Complex64 *alpha, 
	     Complex64 *A, int *lda, Complex64 *x, int *incx, 
	     Complex64 *beta, Complex64 *y, int *incy);

  
  void dsteqr_(char *, int *, double *, double *, double *, int *,
               double *, int *);

  /* Fort Justin */
  void dlartg_(double* F, double* G, double* CS, double* SN, double* R);

  int zgeev_(char *jobvl, 
	    char *jobvr,
	    int  *n,
	    Complex64* A,
	    int* lda,
	    Complex64* w,
	    Complex64* vl,
	    int* ldvl,
	    Complex64* vr,
	    int* ldvr,
	    Complex64* work,
	    int* lwork,
	    double* rwork,
	    int* info);
}

