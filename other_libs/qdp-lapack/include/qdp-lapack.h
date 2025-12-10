// -*- C++ -*-
// $Id: qdp-lapack.h,v 1.8 2009-07-17 14:52:54 bjoo Exp $
/*! \file
 *  \brief QDP interface to Lapack lib
 */

#ifndef __qdp_lapack_h__
#define __qdp_lapack_h__

#include "qdp.h"
using namespace QDP;

namespace QDPLapack
{
  // QDP Interface to Lapack functions
  /*!
   * @{
   */
  using namespace QDP;

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w,
	    multi1d<DComplex>& Work, // Should be length LWork >= max(1,2*N-1)
	    //const int& LWork,
	    multi1d<Double>& RWork // Should be length max(1,3*N-2)
    ) ; 

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w
    ) ;

  int zheev(char& jobz,
	    char& uplo,
	    const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w
    ) ;

  int zgeqrf(const int M, // The std::vector length
	     const int N, // The number of vectors
	     multi2d<DComplex>& A, // the array containing the vectors
	     multi1d<DComplex>& TAU // some strange LAPACK beast
    );

  int zgeev(int n,
	    multi2d<DComplex>& A, 
	    multi1d<DComplex>& evals,
	    multi2d<DComplex>& evecs);

  //! Applies the Q part of a QR factorization 
  /*!
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
   */
  int zunmqr(char& side,//'L' or 'R': apply Q or Q^\dagger from the Left
	     char& trans,//'N' or 'C':  apply Q or Q^\dagger
	     const int M, // The std::vector length
	     const int N,  // The number of vectors
	     multi2d<DComplex>& A, //input MxN matrix
	     multi1d<DComplex>& TAU, // some strange LAPACK beast
	     multi2d<DComplex>& C // input,output
    );    

  int zunmqr2(char& side,//'L' or 'R': apply Q or Q^\dagger from the Left
	     char& trans,//'N' or 'C':  apply Q or Q^\dagger
	     const int M, // The std::vector length
	     const int N,  // The number of vectors
	     const int K, // how many vectors to apply from tau
	     multi2d<DComplex>& A, //input MxN matrix
	     multi1d<DComplex>& TAU, // some strange LAPACK beast
	     multi2d<DComplex>& C // input,output
    );    

  int zunmqrv(char& side,//'L' or 'R': apply Q or Q^\dagger from the Left
	     char& trans,//'N' or 'C':  apply Q or Q^\dagger
	     const int M, // The std::vector length
	     const int K, // how many vectors to apply from tau
	     multi2d<DComplex>& A, //input MxN matrix
	     multi1d<DComplex>& TAU, // some strange LAPACK beast
	     multi1d<DComplex>& C // input,output
    );    

  /* Explicitly form the Q from the Q_R factorization */
  /* Purpose
   *  =======
   *
   *  ZUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
   *  which is defined as the first N columns of a product of K elementary
   *  reflectors of order M
   *
   *        Q  =  H(1) H(2) . . . H(k)
   *
   *  as returned by ZGEQRF.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix Q. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix Q. M >= N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines the
   *          matrix Q. N >= K >= 0.
   *
   *  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
   *          On entry, the i-th column must contain the vector which
   *          defines the elementary reflector H(i), for i = 1,2,...,k, as
   *          returned by ZGEQRF in the first k columns of its array
   *          argument A.
   *          On exit, the M-by-N matrix Q.
   *
   *  LDA     (input) INTEGER
   *          The first dimension of the array A. LDA >= max(1,M).
   *
   *  TAU     (input) COMPLEX*16 array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by ZGEQRF.
   *
   *  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= max(1,N).
   *          For optimum performance LWORK >= N*NB, where NB is the
   *          optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument has an illegal value
   */
  int zungqr(const int M, // The number of rows in Q
	     const int N,  // The number of columns in Q
	     const int K, // The number of reflections in Tau
	     multi2d<DComplex>& A, //input MxN matrix
	     multi1d<DComplex>& TAU // some strange LAPACK beast
	      );    
  
  



  /*--------------------------------------------------------------------
   *  ZHETRF( UPLO, N, A, IPIV)
   *  Double precision Cholesky factorization from Lapack
   *    	A = L*D*L**H
   *
   *  UPLO  'U' stores upper triangular / 'L' stores lower triangular
   *  N     order of A
   *  A     (input)  The matrix 
   *  	    (output) The upper/lower triangular factor stored in A
   *  IPIV  (output) INTEGER array, dimension (N) the pivot array
   *
   *  Work is allocated locally
   *--------------------------------------------------------------------*/
  int zhetrf( char &uplo, const int& n,
	      multi2d<DComplex>& A, 
	      multi1d<int>& ipiv
	      ) ;

/*--------------------------------------------------------------------
   * ZHETRS( UPLO, N, A, IPIV, B)
   *  
   *  Double precision solution of *only NRHS=1* linear system 
   *  using Cholesky factors A = L*D*L**H from ZHETRF
   *
   *  UPLO  'U' stores upper triangular / 'L' stores lower triangular
   *  N     order of A
   *  A     (input) The factor in U or L
   *  IPIV  (input) INTEGER array, dimension (N) the pivot array
   *  B	    one right hand size (NRHS = 1)
   *
   *--------------------------------------------------------------------*/
  int zhetrs( char &uplo, const int& n, 
	      multi2d<DComplex>& A, 
	      multi1d<int>& ipiv,
	      multi1d<DComplex>& B
	      ) ;

  // This is a BLAS routine 
  //     C := alpha*op( A )*op( B ) + beta*C,
  //
  // TRANSA TRANSB 'N' 'T' or 'C' (normal transpose dagger )
  int LatFermMat_x_Mat_cgemm(char& TRANSA,  
			     char& TRANSB,
			     const int& M,  // rows of A
			     const int& N,  // columns of B
			     const int& K,  // columns of A and rows of B
			     const Complex& ALPHA, 
			     const multi1d<LatticeFermion>& A, 
			     // The LDA is known: int *LDA,
			     const multi2d<Complex>& B, const int& LDB,
			     const Complex& BETA,
			     multi1d<LatticeFermion>& C
			     // The LDB is known: int *LDC
			     ) ;

  //Mixed precision
  int LatFermMat_x_Mat_cgemm(char& TRANSA, 
                             char& TRANSB,
                             const int& M, const int& N, const int& K, 
                             const DComplex& ALPHA,
                             const multi1d<LatticeFermion>& A,
                             // The LDA is known: int *LDA,
                             const multi2d<DComplex>& B, const int& LDB, 
                             const Complex& BETA,
                             multi1d<LatticeFermion>& C
                             // The LDB is known: int *LDC
                             ) ;

  /*--------------------------------------------------------------------
   * cgemm     	** this is a BLAS not a Lapack function **
   * 		COMPLEX 
   *            Full interface
   *
   *            C = alpha*A*B + beta*C
   *
   * 		A mxk,	B kxn, 	C mxn 
   * 		(or more accurately op(A),op(B),op(C) have these dimensions)
   * 		
   *--------------------------------------------------------------------*/
  int cgemm(char &transa,
	    char &transb,
	    int m,
	    int n,
	    int k,
	    Complex alpha,
	    multi2d<Complex>& A,  // input
	    int lda,
	    multi2d<Complex>& B,  // input
	    int ldb,
	    Complex beta,
	    multi2d<Complex>& C,  //input/output
	    int ldc   ) ;


/*--------------------------------------------------------------------
   * zgemm     	** this is a BLAS not a Lapack function **
   * 		DOUBLE COMPLEX 
   *            Full interface
   *
   *            C = alpha*A*B + beta*C
   *
   * 		A mxk,	B kxn, 	C mxn
   * 		(or more accurately op(A),op(B),op(C) have these dimensions)
   * 		
   *--------------------------------------------------------------------*/
  int zgemm(char &transa,
	    char &transb,
	    int m,
	    int n,
	    int k,
	    DComplex alpha,
	    multi2d<DComplex>& A,  // input
	    int lda,
	    multi2d<DComplex>& B,  // input
	    int ldb,
	    DComplex beta,
	    multi2d<DComplex>& C,  //input/output
	    int ldc   	    ) ;

  /*--------------------------------------------------------------------
   * CGEMV  	** this is a BLAS not a Lapack function **
   * 		single precision interface for 
   * 			y=A*x
   *  
   *  TRANS   'N' y := A*x,   'T' y := A'*x,   'C' y := conjg( A' )*x 
   *  M     rows of A
   *  N     columns of A
   *  A     The m x n matrix
   *  X     the std::vector of dim (n)
   *  Y     (output) the result
   *
   *--------------------------------------------------------------------*/
  int cgemv(char &trans, 
	    const int& m, const int& n,
	    multi2d<Complex>& A, 
	    multi1d<Complex>& x,
	    multi1d<Complex>& y
	    );

  /*--------------------------------------------------------------------
   * CGEMV  	** this is a BLAS not a Lapack function **
   * 		double precision interface for 
   * 			y=A*x
   *  
   *  TRANS   'N' y := A*x,   'T' y := A'*x,   'C' y := conjg( A' )*x 
   *  M     rows of A
   *  N     columns of A
   *  A     The m x n matrix
   *  X     the std::vector of dim (n)
   *  Y     (output) the result
   *
   *--------------------------------------------------------------------*/
  int zgemv(char &trans, 
	    const int& m, const int& n,
	    multi2d<DComplex>& A, 
	    int lda,
	    multi1d<DComplex>& x,
	    multi1d<DComplex>& y
	    );


  /**
     SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
     *
     *  -- LAPACK routine (version 3.1) --
     *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
     *     November 2006
     *
     *     .. Scalar Arguments ..
     CHARACTER          UPLO
     INTEGER            INFO, LDA, N
     *     ..
     *     .. Array Arguments ..
     COMPLEX*16         A( LDA, * )
     *     ..
     *
     *  Purpose
     *  =======
     *
     *  ZPOTRF computes the Cholesky factorization of a complex Hermitian
     *  positive definite matrix A.
     *
     *  The factorization has the form
     *     A = U**H * U,  if UPLO = 'U', or
     *     A = L  * L**H,  if UPLO = 'L',
     *  where U is an upper triangular matrix and L is lower triangular.
     *
     *  This is the block version of the algorithm, calling Level 3 BLAS.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
     *          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
     *          N-by-N upper triangular part of A contains the upper
     *          triangular part of the matrix A, and the strictly lower
     *          triangular part of A is not referenced.  If UPLO = 'L', the
     *          leading N-by-N lower triangular part of A contains the lower
     *          triangular part of the matrix A, and the strictly upper
     *          triangular part of A is not referenced.
     *
     *          On exit, if INFO = 0, the factor U or L from the Cholesky
     *          factorization A = U**H*U or A = L*L**H.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, the leading minor of order i is not
     *                positive definite, and the factorization could not be
     *                completed.
     *
     **/

  int zpotrf(char &uplo, int N,  multi2d<DComplex>& A, int LDA, int& info) ;
  int zpotrf(char &uplo, multi2d<DComplex>& A, int& info) ;
  
  int cpotrf(char &uplo, int N,  multi2d<Complex>& A, int LDA, int& info) ;
  int cpotrf(char &uplo, multi2d<Complex>& A, int& info) ;

  int zpotrf(char &uplo, int N,  multi1d<DComplex>& A, int LDA, int& info) ;
  
  int cpotrf(char &uplo, int N,  multi1d<Complex>& A, int LDA, int& info) ;

  /**
     SUBROUTINE ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
     *
     *  -- LAPACK routine (version 3.1) --
     *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
     *     November 2006
     *
     *     .. Scalar Arguments ..
     CHARACTER          UPLO
     INTEGER            INFO, LDA, LDB, N, NRHS
     *     ..
     *     .. Array Arguments ..
     COMPLEX*16         A( LDA, * ), B( LDB, * )
     *     ..
     *
     *  Purpose
     *  =======
     *
     *  ZPOTRS solves a system of linear equations A*X = B with a Hermitian
     *  positive definite matrix A using the Cholesky factorization
     *  A = U**H*U or A = L*L**H computed by ZPOTRF.
     *
     *  Arguments
     *  =========
     *
     *  UPLO    (input) CHARACTER*1
     *          = 'U':  Upper triangle of A is stored;
     *          = 'L':  Lower triangle of A is stored.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A.  N >= 0.
     *
     *  NRHS    (input) INTEGER
     *          The number of right hand sides, i.e., the number of columns
     *          of the matrix B.  NRHS >= 0.
     *
     *  A       (input) COMPLEX*16 array, dimension (LDA,N)
     *          The triangular factor U or L from the Cholesky factorization
     *          A = U**H*U or A = L*L**H, as computed by ZPOTRF.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
     *          On entry, the right hand side matrix B.
     *          On exit, the solution matrix X.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of the array B.  LDB >= max(1,N).
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *
     *  =====================================================================
     **/

  int zpotrs(char &uplo, int N, int nrhs,  
	     multi2d<DComplex>& A, int LDA, 
	     multi2d<DComplex>& B, int LDB, int& info) ;

  int zpotrs(char &uplo, multi2d<DComplex>& A,  multi2d<DComplex>& B,  
	     int& info) ;

  int cpotrs(char &uplo, int N, int nrhs,  
	     multi2d<Complex>& A, int LDA, 
	     multi2d<Complex>& B, int LDB, int& info) ;

  int cpotrs(char &uplo, multi2d<Complex>& A,  multi2d<Complex>& B,  
	     int& info) ;


  int zpotrs(char &uplo, int N, int nrhs,  
	     multi1d<DComplex>& A, int LDA, 
	     multi2d<DComplex>& B, int LDB, int& info) ;

  int cpotrs(char &uplo, int N, int nrhs,  
	     multi1d<Complex>& A, int LDA, 
	     multi2d<Complex>& B, int LDB, int& info) ;

  
  void dsteqr(char *, int *, double *, double *, double *, int *,
              double *, int *);


  //for Justin
  void dlartg(double& F, double& G, double& CS, double& SN, double& R);

  // LU Factorize a GENERAL DComplex Matrix
  //  M is the number of rows
  //  N is the number of columns
  //  LDA is the leading dimension of A
  //  IPIV is the Integer array of pivots
  //  Info is a return value
  int zgetrf(int M, int N, multi2d<Complex64>& A, int LDA,
	     multi1d<int>& ipiv, int& info);

  void zgetrs(char& trans, int n, int rhs,
	     multi2d<Complex64>& A, int lda, 
	     multi1d<int>& ipiv, 
	     multi1d<Complex64>& B,
	     int ldb, int& info);

} // namespace QDPLapack

#endif 
