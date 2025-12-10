// -*- C++ -*-
// $Id: fortran_lapack.cc,v 1.9 2009-07-31 14:09:31 bjoo Exp $
/*! \file
 *  \brief QDP interface to Lapack lib using c-lapack
 */

#include "qdp-lapack.h"
#include "fortran_lapack.h"

namespace QDPLapack
{
  // hide lapack's temporaries...
  int zheev(char& jobz,
	    char& uplo,
	    const int& N, // These are the dimensions of the array A
	    multi2d<Complex64>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w)
  {
    /**
       char c_jobz = *jobz;
       char c_uplo = *uplo;
    **/

    int LWork = 2*N-1;
    multi1d<Complex64> Work(LWork);
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
      QDPIO::cerr<<"Lapack::zheev returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;
  }

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<Complex64>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w,
	    multi1d<Complex64>& Work, // Should be length LWork >= max(1,2*N-1)
	    //const int& LWork,
	    multi1d<Double>& RWork // Should be length max(1,3*N-2)
    )
  {
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
      QDPIO::cerr<<"Lapack::zheev returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;

  }

  // hide lapack's temporaries...
  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<Complex64>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w)
  {
    int N = A.size1();
    int LWork = 2*N-1;
    multi1d<Complex64> Work(LWork);
    multi1d<Double> RWork(3*N-2);
    
    return zheev(jobz,uplo,A,w,Work,RWork) ;
  }

  int zgeqrf(const int M, // The std::vector length
	     const int N, // The number of vectors
	     multi2d<Complex64>& A, // the array containing the vectors
	     multi1d<Complex64>& TAU // some strange LAPACK beast
    )
  {
    if(N>M)
      TAU.resize(M);
    else
      TAU.resize(N);
	
    int LWork = N ;
    multi1d<Complex64> Work(LWork);
    

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
      QDPIO::cerr<<"Lapack::zgeqrf returned with exit code: "<<info<<std::endl ;
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
	     multi2d<Complex64>& A, //input
	     multi1d<Complex64>& TAU, // some strange LAPACK beast
	     multi2d<Complex64>& C //input/output
    )
  {
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
       std::cout<<"ZUNMQR->side  : "<<side<<std::endl  ;
       std::cout<<"ZUNMQR->trans : "<<trans<<std::endl ;
       std::cout<<"ZUNMQR->M     : "<<M<<std::endl ;
       std::cout<<"ZUNMQR->N     : "<<N<<std::endl ;
       std::cout<<"ZUNMQR->K     : "<<K<<std::endl ;
       std::cout<<"ZUNMQR->lda   : "<<lda<<std::endl ;
       std::cout<<"ZUNMQR->ldc   : "<<ldc<<std::endl ;
       std::cout<<"ZUNMQR->LWork : "<<LWork<<std::endl ;
    **/

    multi1d<Complex64> Work(LWork);
    
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
      QDPIO::cerr<<"Lapack::zunmqr returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;
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
    /**
       char c_side = *side ;
       char c_trans = *trans ;
    **/

    int lda = A.size1();
    int ldc = C.size1();
    int LWork = -1 ;
    if(side == 'R')
      LWork = M ;
    if(side == 'L')
      LWork = N ;

    /**
       std::cout<<"ZUNMQR->side  : "<<side<<std::endl  ;
       std::cout<<"ZUNMQR->trans : "<<trans<<std::endl ;
       std::cout<<"ZUNMQR->M     : "<<M<<std::endl ;
       std::cout<<"ZUNMQR->N     : "<<N<<std::endl ;
       std::cout<<"ZUNMQR->K     : "<<K<<std::endl ;
       std::cout<<"ZUNMQR->lda   : "<<lda<<std::endl ;
       std::cout<<"ZUNMQR->ldc   : "<<ldc<<std::endl ;
       std::cout<<"ZUNMQR->LWork : "<<LWork<<std::endl ;
    **/

    multi1d<Complex64> Work(LWork);
    
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
      QDPIO::cerr<<"Lapack::zunmqr returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;
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
    /**
       char c_side = *side ;
       char c_trans = *trans ;
    **/

    int lda = A.size1();
    int ldc = C.size();
    int LWork = -1 ;
    int N=1;
    if(side == 'R')
      LWork = M ;
    if(side == 'L')
      LWork = N ;

    /**
       std::cout<<"ZUNMQR->side  : "<<side<<std::endl  ;
       std::cout<<"ZUNMQR->trans : "<<trans<<std::endl ;
       std::cout<<"ZUNMQR->M     : "<<M<<std::endl ;
       std::cout<<"ZUNMQR->N     : "<<N<<std::endl ;
       std::cout<<"ZUNMQR->K     : "<<K<<std::endl ;
       std::cout<<"ZUNMQR->lda   : "<<lda<<std::endl ;
       std::cout<<"ZUNMQR->ldc   : "<<ldc<<std::endl ;
       std::cout<<"ZUNMQR->LWork : "<<LWork<<std::endl ;
    **/

    multi1d<Complex64> Work(LWork);
    
    int info ;
    int r = zunmqr_(&side, &trans,
		    (int *)&M, 
		    (int *)&N, 
		    (int *)&K, 
		    &A(0,0),
		    &lda, 
		    &TAU[0], 
		    &C[0],
		    &ldc,
		    &Work[0],
		    &LWork,
		    &info);

    if(info){
      QDPIO::cerr<<"QDPLapack::zunmqrv returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;
  }
    
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
	     )
  {

    int lda = A.size1();
    int info;
    multi1d<Complex64> work(1); 
    int lwork=-1;
    
    // Query workspace
    zungqr_(&M, &N, &K, &A(0,0), &lda, &TAU[0], &work[0],&lwork, &info);
    if( info != 0) { 
      QDPIO::cout << "Lapack failed; Info: argument " << (-info) << " has illegal value " <<std::endl;
      QDP_abort(1);
    }

    lwork=toInt(real(work[0]));
    work.resize(lwork);
    // Do the work
    zungqr_(&M, &N, &K, &A(0,0), &lda, &TAU[0], &work[0],&lwork, &info);
    if( info != 0) { 
      QDPIO::cout << "Lapack failed; Info: argument " << (-info) << " has illegal value " <<std::endl;
      QDP_abort(1);
    }

  }   

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
	      multi2d<Complex64>& A, 
	      multi1d<int>& ipiv
	      )
  {
    int info;
    int lda = A.size1() ; 
    int LWork = n;
    multi1d<Complex64> Work(LWork);
    
    int r = zhetrf_( &uplo, (int *)&n, &A(0,0), &lda, 
		     &ipiv[0],  &Work[0], &LWork, &info );
    
    if(info){
      QDPIO::cerr<<"Lapack::zhetrf returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;
  }

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
	      multi2d<Complex64>& A, 
	      multi1d<int>& ipiv,
	      multi1d<Complex64>& B
	    )
  {
    int NRHS(1);
    int info;
    int lda = A.size1() ; 
    int ldb = B.size() ; 

    int r = zhetrs_( &uplo, (int *)&n, &NRHS,
		     &A(0,0), &lda, 
		     &ipiv[0], 
		     &B[0], &ldb,
		     &info );
    
    if(info){
      QDPIO::cerr<<"Lapack::zhetrs returned with exit code: "<<info<<std::endl ;
      exit(1) ;
    }

    return r ;
  }


  // Interfaces to BLAS start here
  // Explicit Single Precision....
  int LatFermMat_x_Mat_cgemm(char& TRANSA, 
                             char& TRANSB,
                             const int& M, const int& N, const int& K, 
                             const Complex32& ALPHA,
                             const multi1d<LatticeFermionF3>& A,
                             // The LDA is known: int *LDA,
                             const multi2d<Complex32>& B, const int& LDB, 
                             const Complex32& BETA,
                             multi1d<LatticeFermionF3>& C
                             // The LDB is known: int *LDC
                             ){

    int LDA,LDC ;

    if(N==1)
      LDA = Layout::sitesOnNode() * Nc *Ns ;
    else
      LDA = (Complex32 *)&A[1].elem(0).elem(0).elem(0) - 
            (Complex32 *)&A[0].elem(0).elem(0).elem(0) ;

    C.resize(N) ;
    if(N==1)
      LDC = Layout::sitesOnNode() * Nc *Ns ;
    else
      LDC = (Complex32 *)&C[1].elem(0).elem(0).elem(0) - 
	    (Complex32 *)&C[0].elem(0).elem(0).elem(0) ;

    //DOES C NEED TO BE ZEROED out ?
    //for(int i(0);i<N;i++) C[i] = zero ;

    /**
       QDPIO::cout<<"Leading dim is LDA= "<<LDA<<std::endl ;
       int ttLD = Layout::sitesOnNode() * Nc *Ns ;
       QDPIO::cout<<"  LatticeFerion Data size= "<<ttLD<<std::endl ;
       QDPIO::cout<<"   extra stuff= "<<LDA-ttLD<<std::endl ;
       QDPIO::cout<<"Leading dim is LDC= "<<LDC<<std::endl ;
       QDPIO::cout<<"   extra stuff= "<<LDC-ttLD<<std::endl ;
    **/

    return cgemm_(&TRANSA,&TRANSB,(int *)&M,(int *)&N,(int *)&K,
		 (Complex32 *)&ALPHA,
		 (Complex32 *)&A[0].elem(0).elem(0).elem(0),&LDA,
		 (Complex32 *)&B(0,0),(int *)&LDB,
		 (Complex32 *)&BETA,
		 (Complex32 *)&C[0].elem(0).elem(0).elem(0),&LDC);
  }

  // Interfaces to BLAS start here
  // Explicit Double Precision....
  int LatFermMat_x_Mat_cgemm(char& TRANSA, 
                             char& TRANSB,
                             const int& M, const int& N, const int& K, 
                             const Complex64& ALPHA,
                             const multi1d<LatticeFermionD3>& A,
                             // The LDA is known: int *LDA,
                             const multi2d<Complex64>& B, const int& LDB, 
                             const Complex64& BETA,
                             multi1d<LatticeFermionD3>& C
                             // The LDB is known: int *LDC
                             ){

    int LDA,LDC ;

    if(N==1)
      LDA = Layout::sitesOnNode() * Nc *Ns ;
    else
      LDA = (Complex64 *)&A[1].elem(0).elem(0).elem(0) - 
            (Complex64 *)&A[0].elem(0).elem(0).elem(0) ;

    C.resize(N) ;
    if(N==1)
      LDC = Layout::sitesOnNode() * Nc *Ns ;
    else
      LDC = (Complex64 *)&C[1].elem(0).elem(0).elem(0) - 
	    (Complex64 *)&C[0].elem(0).elem(0).elem(0) ;

    //DOES C NEED TO BE ZEROED out ?
    //for(int i(0);i<N;i++) C[i] = zero ;

    /**
       QDPIO::cout<<"Leading dim is LDA= "<<LDA<<std::endl ;
       int ttLD = Layout::sitesOnNode() * Nc *Ns ;
       QDPIO::cout<<"  LatticeFerion Data size= "<<ttLD<<std::endl ;
       QDPIO::cout<<"   extra stuff= "<<LDA-ttLD<<std::endl ;
       QDPIO::cout<<"Leading dim is LDC= "<<LDC<<std::endl ;
       QDPIO::cout<<"   extra stuff= "<<LDC-ttLD<<std::endl ;
    **/

    return zgemm_(&TRANSA,&TRANSB,(int *)&M,(int *)&N,(int *)&K,
		 (Complex64 *)&ALPHA,
		 (Complex64 *)&A[0].elem(0).elem(0).elem(0),&LDA,
		 (Complex64 *)&B(0,0),(int *)&LDB,
		 (Complex64 *)&BETA,
		 (Complex64 *)&C[0].elem(0).elem(0).elem(0),&LDC);
  }


#if 0  
  // I don't know if this was really needed ever, or whether
  // It was an attempt to make a double precision build work.
  // I am commenting it out now...
  
  //Mixed precision
  int LatFermMat_x_Mat_cgemm(char& TRANSA, 
                             char& TRANSB,
                             const int& M, const int& N, const int& K, 
                             const Complex64& ALPHA,
                             const multi1d<LatticeFermion>& A,
                             // The LDA is known: int *LDA,
                             const multi2d<Complex64>& B, const int& LDB, 
                             const Complex64& BETA,
                             multi1d<LatticeFermion>& C
                             // The LDB is known: int *LDC
                             ){
    Complex alpha = ALPHA ;
    Complex beta = BETA ;
    
    multi2d<Complex> b(B.size2(),B.size1());
    for(int i(0);i<B.size2();i++)
      for(int j(0);j<B.size1();j++)
	b(i,j) = B(i,j);
    
    return LatFermMat_x_Mat_cgemm(TRANSA, TRANSB,M, N, K, alpha, A, b, LDB, 
                             beta, C  );
  }
#endif


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
	    Complex32 alpha,
	    multi2d<Complex32>& A,  // input
	    int lda,
	    multi2d<Complex32>& B,  // input
	    int ldb,
	    Complex32 beta,
	    multi2d<Complex32>& C,  //input/output
	    int ldc   
    )
  {
    int r = cgemm_(&transa,&transb, &m, &n, &k, &alpha,
		   &A(0,0),&lda, &B(0,0),&ldb, &beta, &C(0,0),&ldc
		   );
  }

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
    int r = zgemm_(&transa,&transb, &m, &n, &k, &alpha,
		   &A(0,0),&lda, &B(0,0),&ldb, &beta, &C(0,0),&ldc
		  );
  }

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
	    multi2d<Complex32>& A, 
	    multi1d<Complex32>& x,
	    multi1d<Complex32>& y
	    )
  {
    int incx(1), incy(1);
    Complex32 alpha(1.0), beta(0.0);
    int lda = A.size1() ; 

    return cgemv_(&trans, (int *) &m, (int *) &n, &alpha, 
		    (Complex32 *) &A(0,0), &lda, 
		    &x[0], &incx, &beta, 
		    &y[0], &incy);
  }


  /*--------------------------------------------------------------------
   * ZGEMV  	** this is a BLAS not a Lapack function **
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
	    multi2d<Complex64>& A, 
	    int lda,
	    multi1d<Complex64>& x,
	    multi1d<Complex64>& y
	    )
  {
    int incx(1), incy(1);
    Complex64 alpha(1.0), beta(0.0);
    return zgemv_(&trans, (int *) &m, (int *) &n, &alpha, 
		    (Complex64 *) &A(0,0), &lda, 
		    &x[0], &incx, &beta, 
		    &y[0], &incy);
  }


  int zpotrf(char &uplo, int N,  multi2d<Complex64>& A, int LDA, int& info){
    
    return zpotrf_(&uplo, &N, 
		   (Complex64 *) &A(0,0), &LDA, &info);

  }
  int zpotrf(char &uplo, multi2d<Complex64>& A, int& info){
    int LDA = A.size1() ; // Assumes square matrix
    
    return zpotrf(uplo,LDA,A,LDA,info);
    
  }


  int cpotrf(char &uplo, int N,  multi2d<Complex32>& A, int LDA, int& info){
    
    return cpotrf_(&uplo, &N, 
		   (Complex32 *) &A(0,0), &LDA, &info);

  }
  int cpotrf(char &uplo, multi2d<Complex32>& A, int& info){
    int LDA = A.size1() ; // Assumes square matrix
    
    return cpotrf(uplo,LDA,A,LDA,info);
    
  }



  int zpotrf(char &uplo, int N,  multi1d<Complex64>& A, int LDA, int& info){
    
    return zpotrf_(&uplo, &N, (Complex64 *) &A[0], &LDA, &info);
  }


  int cpotrf(char &uplo, int N,  multi1d<Complex32>& A, int LDA, int& info){
    
    return cpotrf_(&uplo, &N, (Complex32 *) &A[0], &LDA, &info);
  }

  
  int zpotrs(char &uplo, int N, int nrhs,  
	     multi2d<Complex64>& A, int LDA, 
	     multi2d<Complex64>& B, int LDB, int& info){

    return zpotrs_(&uplo, &N, &nrhs, 
		   (Complex64 *)&A(0,0), &LDA,
		   (Complex64 *)&B(0,0), &LDB,
		   &info);
  }
  
  int zpotrs(char &uplo, multi2d<Complex64>& A,  multi2d<Complex64>& B,  
	     int& info){
    
    int LDA = A.size1() ; // Assumes square matrix
    int LDB = B.size1() ; // second index is fast and its size is size1
    int nrhs = B.size2() ;
    
    zpotrs(uplo,LDA,nrhs,A,LDA,B,LDB,info);
    
  }

  

  int cpotrs(char &uplo, int N, int nrhs,  
	     multi2d<Complex32>& A, int LDA, 
	     multi2d<Complex32>& B, int LDB, int& info){

    return cpotrs_(&uplo, &N, &nrhs, 
		   (Complex32 *)&A(0,0), &LDA,
		   (Complex32 *)&B(0,0), &LDB,
		   &info);
  }
  
  int cpotrs(char &uplo, multi2d<Complex32>& A,  multi2d<Complex32>& B,  
	     int& info){
    
    int LDA = A.size1() ; // Assumes square matrix
    int LDB = B.size1() ; // second index is fast and its size is size1
    int nrhs = B.size2() ;
    
    cpotrs(uplo,LDA,nrhs,A,LDA,B,LDB,info);
    
  }
  

  int zpotrs(char &uplo, int N, int nrhs,  
	     multi1d<Complex64>& A, int LDA, 
	     multi2d<Complex64>& B, int LDB, int& info){
    
    return zpotrs_(&uplo, &N, &nrhs, 
		   (Complex64 *)&A[0], &LDA,
		   (Complex64 *)&B(0,0), &LDB,
		   &info);
  }
  
  
  int cpotrs(char &uplo, int N, int nrhs,  
	     multi1d<Complex32>& A, int LDA, 
	     multi2d<Complex32>& B, int LDB, int& info){

    return cpotrs_(&uplo, &N, &nrhs, 
		   (Complex32 *)&A[0], &LDA,
		   (Complex32 *)&B(0,0), &LDB,
		   &info);
  }


  void dsteqr(char *a, int *b, double *c, double *d, double *e, int *f,
              double *g, int *h)
  {
    dsteqr_(a,b,c,d,e,f,g,h);
  }


  // Pass by reference, so we can take addresses of the reference args.
  // for the fortran cross call.
  void dlartg(double& F, double& G, double& CS, double& SN, double& R)
  {
    dlartg_(&F,&G,&CS,&SN,&R);
  }


  // ZGETRF
  int zgetrf(int M, int N, multi2d<Complex64>& A, int LDA,
	      multi1d<int>& ipiv, int& info)
  {
    int r = zgetrf_(&M, &N, (Complex64*)&A(0,0), &LDA, &ipiv[0], &info);

    if( info != 0 ) { 
      if( info < 0 ) {
	QDPIO::cout << "ZGETRF: argument"<< (-info) << " has illegal value" << std::endl;
	QDP_abort(1);
      }
      if( info > 0 ) { 
	QDPIO::cout << "ZGETRF: U(" << info <<","<< info <<") is exactly zero making U singular" << std::endl;
	QDP_abort(1);
      }
    }

    return r;
  }
  // ZGETRS
  
  void zgetrs(char &trans, int n, int nrhs,
	     multi2d<Complex64>& A, int lda, 
	     multi1d<int>& ipiv,  
	     multi1d<Complex64>& B,
	     int ldb, int& info) 
  {
    zgetrs_(&trans,
	    &n,
	    &nrhs,
	    &A(0,0),
	    &lda,
	    &ipiv[0],
	    &B[0],
	    &ldb,
	    &info);
    
    if (info) {
      QDPIO::cerr << "qdp_lapack: Linear system solve by ZGETRS failed.\n" << std::endl;
    }
  }  

  int zgeev(int n,
	    multi2d<DComplex>& A, 
	    multi1d<DComplex>& evals,
	    multi2d<DComplex>& evecs)
  {
    char jobvl='N';
    char jobvr='V';
    int lda = A.size1();
    int ldvr= n;

    int lwork=4*n;
    multi2d<Complex64> vl(n,n);
  
    multi1d<Complex64> work(lwork); // Work array
    multi1d<double> rwork(2*n);
    int info;

    zgeev_(&jobvl,
	   &jobvr,
	   &n,
	   &A(0,0),
	   &lda,
	   &evals[0],
	   &vl(0,0), 
	   &n, 
	   &evecs(0,0),
	   &ldvr,
	   &work[0],
	   &lwork,
	   &rwork[0],
	   &info);

    if ( info ) { 
      if( info < 0 ) { 
	QDPIO::cerr << "QDPLapack::zgeev: argument "<< (-info) << " to LAPACK ZGEEV had illegal value" << std::endl;
      }
      else { 
	QDPIO::cerr << "QDPLapack::zgeev: LAPACK ZGEEV Eigensolver failed" << std::endl;
      }
      QDP_abort(1);
    }


      

  }

  
} // namespace QDPLapack


