/*******************************************************************************
 *   Adapted from PRIMME Feb 7, 2008
 *
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
 *
 ******************************************************************************/

#ifndef NUMERICAL_H
#define NUMERICAL_H

#if !defined(NUM_SUN) && !defined(NUM_IBM) && !defined(NUM_CRAY)
#define NUM_SUN
#endif

#ifdef NUM_SUN


#define ZDOTCSUB  zdotcsub_
#define BLAS_ZCOPY  zcopy_
#define BLAS_ZSWAP  zswap_
#define BLAS_ZGEMM  zgemm_
#define BLAS_ZHEMM  zhemm_
#define BLAS_ZAXPY  zaxpy_
#define BLAS_ZGEMV  zgemv_
#define BLAS_ZSCAL  zscal_
#define BLAS_ZLARNV zlarnv_
#define BLAS_ZHEEV  zheev_
#define BLAS_ZHETRF zhetrf_
#define BLAS_ZHETRS zhetrs_
#define BLAS_ZGEQRF zgeqrf_
#define BLAS_ZUNMQR zunmqr_


#define CDOTCSUB  cdotcsub_
#define BLAS_CCOPY  ccopy_
#define BLAS_CSWAP  cswap_
#define BLAS_CGEMM  cgemm_
#define BLAS_CHEMM  chemm_
#define BLAS_CAXPY  caxpy_
#define BLAS_CGEMV  cgemv_
#define BLAS_CSCAL  cscal_
#define BLAS_CLARNV clarnv_
#define BLAS_CHEEV  cheev_
#define BLAS_CHETRF chetrf_
#define BLAS_CHETRS chetrs_
#define BLAS_CGEQRF cgeqrf_
#define BLAS_CUNMQR cunmqr_
#define BLAS_CPOTRF cpotrf_
#define BLAS_CPOTRS cpotrs_

#define BLAS_DCOPY  dcopy_
#define BLAS_DSWAP  dswap_
#define BLAS_DGEMM  dgemm_
#define BLAS_DSYMM  dsymm_
#define BLAS_DAXPY  daxpy_
#define BLAS_DGEMV  dgemv_
#define BLAS_DDOT   ddot_
#define BLAS_DSCAL  dscal_
#define BLAS_DLARNV dlarnv_
#define BLAS_DSYEV  dsyev_
#define BLAS_DSYTRF dsytrf_
#define BLAS_DSYTRS dsytrs_

#define BLAS_SCOPY  scopy_
#define BLAS_SSWAP  sswap_
#define BLAS_SGEMM  sgemm_
#define BLAS_SSYMM  ssymm_
#define BLAS_SAXPY  saxpy_
#define BLAS_SGEMV  sgemv_
#define BLAS_SDOT   sdot_
#define BLAS_SSCAL  sscal_
#define BLAS_SLARNV slarnv_
#define BLAS_SSYEV  ssyev_
#define BLAS_SSYTRF ssytrf_
#define BLAS_SSYTRS ssytrs_

/*Functions added for use with EIGBICG*/
/**********************************/
#define BLAS_CGEEV  cgeev_
#define BLAS_ZGEEV  zgeev_
#define BLAS_CGESV  cgesv_
#define BLAS_ZGESV  zgesv_
#define BLAS_CGETRF cgetrf_
#define BLAS_CGETRS cgetrs_
#define BLAS_ZGETRF zgetrf_
#define BLAS_ZGETRS zgetrs_
/*********************************/



#elif defined(NUM_IBM)



#define ZDOTCSUB  zdotcsub
#define BLAS_ZCOPY  zcopy
#define BLAS_ZSWAP  zswap
#define BLAS_ZGEMM  zgemm
#define BLAS_ZHEMM  zhemm
#define BLAS_ZAXPY  zaxpy
#define BLAS_ZGEMV  zgemv
#define BLAS_ZSCAL  zscal
#define BLAS_ZLARNV zlarnv
#define BLAS_ZHEEV  zheev
#define BLAS_ZHETRF zhetrf
#define BLAS_ZHETRS zhetrs
#define BLAS_ZGEQRF zgeqrf
#define BLAS_ZUNMQR zunmqr

#define CDOTCSUB  cdotcsub
#define BLAS_CCOPY  ccopy
#define BLAS_CSWAP  cswap
#define BLAS_CGEMM  cgemm
#define BLAS_CHEMM  chemm
#define BLAS_CAXPY  caxpy
#define BLAS_CGEMV  cgemv
#define BLAS_CSCAL  cscal
#define BLAS_CLARNV clarnv
#define BLAS_CHEEV  cheev
#define BLAS_CHETRF chetrf
#define BLAS_CHETRS chetrs
#define BLAS_CGEQRF cgeqrf
#define BLAS_CUNMQR cunmqr

#define BLAS_DCOPY  dcopy
#define BLAS_DSWAP  dswap
#define BLAS_DGEMM  dgemm
#define BLAS_DSYMM  dsymm
#define BLAS_DAXPY  daxpy
#define BLAS_DGEMV  dgemv
#define BLAS_DDOT   ddot
#define BLAS_DSCAL  dscal
#define BLAS_DLARNV dlarnv
#define BLAS_DSYEV  dsyev
#define BLAS_DSYTRF dsytrf
#define BLAS_DSYTRS dsytrs


#define BLAS_SCOPY  scopy
#define BLAS_SSWAP  sswap
#define BLAS_SGEMM  sgemm
#define BLAS_SSYMM  ssymm
#define BLAS_SAXPY  saxpy
#define BLAS_SGEMV  sgemv
#define BLAS_SDOT   sdot
#define BLAS_SSCAL  sscal
#define BLAS_SLARNV slarnv
#define BLAS_SSYEV  ssyev
#define BLAS_SSYTRF ssytrf
#define BLAS_SSYTRS ssytrs

/*Functions added for use with EIGBICG*/
/**********************************/
#define BLAS_CGEEV  cgeev
#define BLAS_ZGEEV  zgeev
#define BLAS_CGESV  cgesv
#define BLAS_ZGESV  zgesv
#define BLAS_CGETRF cgetrf
#define BLAS_CGETRS cgetrs
#define BLAS_ZGETRF zgetrf
#define BLAS_ZGETRS zgetrs
/*********************************/



#ifdef NUM_ESSL
#include <essl.h>
#endif

#elif defined(NUM_CRAY)
#include <fortran.h>
#include <string.h>



#define ZDOTCSUB  zdotcsub
#define BLAS_ZCOPY  zcopy
#define BLAS_ZSWAP  zswap
#define BLAS_ZGEMM  zgemm
#define BLAS_ZHEMM  zhemm
#define BLAS_ZAXPY  zaxpy
#define BLAS_ZGEMV  zgemv
#define BLAS_ZSCAL  zscal
#define BLAS_ZLARNV zlarnv
#define BLAS_ZHEEV  zheev
#define BLAS_ZHETRF zhetrf
#define BLAS_ZHETRS zhetrs
#define BLAS_ZGEQRF zgeqrf
#define BLAS_ZUNMQR zunmqr

#define BLAS_DCOPY  SCOPY
#define BLAS_DSWAP  SSWAP
#define BLAS_DGEMM  SGEMM
#define BLAS_DSYMM  DSYMM
#define BLAS_DAXPY  SAXPY
#define BLAS_DGEMV  SGEMV
#define BLAS_DDOT   SDOT
#define BLAS_DLAMCH SLAMCH
#define BLAS_DSCAL  SSCAL
#define BLAS_DLARNV SLARNV
#define BLAS_DSYEV  SSYEV
#define BLAS_DSYTRF DSYTRF
#define BLAS_DSYTRS DSYTRS

/*Functions added for use with EIGBICG*/
/**********************************/
#define BLAS_CGEEV  cgeev
#define BLAS_ZGEEV  zgeev
#define BLAS_CGESV  cgesv
#define BLAS_ZGESV  zgesv
#define BLAS_CGETRF cgetrf
#define BLAS_CGETRS cgetrs
#define BLAS_ZGETRF zgetrf
#define BLAS_ZGETRS zgetrs
/*********************************/



#endif
#ifdef Cplusplus
extern "C" {
#endif /* Cplusplus */

#ifndef NUM_CRAY
	/* If not CRAY then use the following headers */

void BLAS_DCOPY(int *n, double *x, int *incx, double *y, int *incy);
void BLAS_DSWAP(int *n, double *x, int *incx, double *y, int *incy);
void BLAS_DGEMM(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void BLAS_DSYMM(char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void BLAS_DAXPY(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void BLAS_DGEMV(char *transa, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
double BLAS_DDOT(int *n, double *x, int *incx, double *y, int *incy);
double BLAS_DLAMCH(char *cmach);
void BLAS_DSCAL(int *n, double *alpha, double *x, int *incx);
void BLAS_DLARNV(int *idist, int *iseed, int *n, double *x);
void BLAS_DSYEV(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *ldwork, int *info);
void BLAS_DSYTRF(char *uplo, int *n, double *a, int *lda, int *ipivot, double *work, int *ldwork, int *info);
void BLAS_DSYTRS(char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipivot, double *b, int *ldb, int *info);

void BLAS_SCOPY(int *n, float *x, int *incx, float *y, int *incy);
void BLAS_SSWAP(int *n, float *x, int *incx, float *y, int *incy);
void BLAS_SGEMM(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
void BLAS_SSYMM(char *side, char *uplo, int *m, int *n, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
void BLAS_SAXPY(int *n, float *alpha, float *x, int *incx, float *y, int *incy);
void BLAS_SGEMV(char *transa, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
float BLAS_SDOT(int *n, float *x, int *incx, float *y, int *incy);
float BLAS_SLAMCH(char *cmach);
void BLAS_SSCAL(int *n, float *alpha, float *x, int *incx);
void BLAS_SLARNV(int *idist, int *iseed, int *n, float *x);
void BLAS_SSYEV(char *jobz, char *uplo, int *n, float *a, int *lda, float *w, float *work, int *ldwork, int *info);
void BLAS_SSYTRF(char *uplo, int *n, float *a, int *lda, int *ipivot, float *work, int *ldwork, int *info);
void BLAS_SSYTRS(char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipivot, float *b, int *ldb, int *info);

void   BLAS_ZCOPY(int *n, void *x, int *incx, void *y, int *incy);
void   BLAS_ZSWAP(int *n, void *x, int *incx, void *y, int *incy);
void   BLAS_ZGEMM(char *transa, char *transb, int *m, int *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   BLAS_ZHEMM(char *side, char *uplo, int *m, int *n, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   BLAS_ZAXPY(int *n, void *alpha, void *x, int *incx, void *y, int *incy);
void   BLAS_ZGEMV(char *transa, int *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
void   BLAS_ZSCAL(int *n, void *alpha, void *x, int *incx);
void   BLAS_ZLARNV(int *idist, int *iseed, int *n, void *x);
void   BLAS_ZHEEV(char *jobz, char *uplo, int *n, void *a, int *lda, double *w, void *work, int *ldwork, double *rwork, int *info);
void   BLAS_ZHETRF(char *uplo, int *n, void *a, int *lda, int *ipivot, void *work, int *ldwork, int *info);
void   BLAS_ZHETRS(char *uplo, int *n, int *nrhs, void *a, int *lda, int *ipivot, void *b, int *ldb, int *info);
void   BLAS_ZGEQRF(int *, int *, void *, int *, void *, void *, int *, int *);
void   BLAS_ZUNMQR(char *,char *,int *,int *,int *,void *,int *, void *, void *, int *, void *,int *, int * );
void   ZDOTCSUB(void *dot, int *n, void *x, int *incx, void *y, int *incy);

void   BLAS_CCOPY(int *n, void *x, int *incx, void *y, int *incy);
void   BLAS_CSWAP(int *n, void *x, int *incx, void *y, int *incy);
void   BLAS_CGEMM(char *transa, char *transb, int *m, int *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   BLAS_CHEMM(char *side, char *uplo, int *m, int *n, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   BLAS_CAXPY(int *n, void *alpha, void *x, int *incx, void *y, int *incy);
void   BLAS_CGEMV(char *transa, int *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
void   BLAS_CSCAL(int *n, void *alpha, void *x, int *incx);
void   BLAS_CLARNV(int *idist, int *iseed, int *n, void *x);
void   BLAS_CHEEV(char *jobz, char *uplo, int *n, void *a, int *lda, float *w, void *work, int *ldwork, float *rwork, int *info);
void   BLAS_CHETRF(char *uplo, int *n, void *a, int *lda, int *ipivot, void *work, int *ldwork, int *info);
void   BLAS_CHETRS(char *uplo, int *n, int *nrhs, void *a, int *lda, int *ipivot, void *b, int *ldb, int *info);
void   BLAS_CGEQRF(int *, int *, void *, int *, void *, void *, int *, int *);
void   BLAS_CUNMQR(char *,char *,int *,int *,int *,void *,int *, void *, void *, int *, void *,int *, int * );
void   BLAS_CHETRF(char *uplo, int *n, void *a, int *lda, int *ipivot, void *work, int *ldwork, int *info);
void   BLAS_CHETRS(char *uplo, int *n, int *nrhs, void *a, int *lda, int *ipivot, void *b, int *ldb, int *info);
void   BLAS_CUNMQR(char *,char *,int *,int *,int *,void *,int *, void *, void *, int *, void *,int *, int * );
void   BLAS_CPOTRF(char *uplo, int *N, void *A, int *LDA, int *info );
void   BLAS_CPOTRS(char *uplo, int *N, int *NRHS, Complex_C *A, int *LDA, Complex_C *B,int *LDB,int *INFO);
void   CDOTCSUB(void *dot, int *n, void *x, int *incx, void *y, int *incy);


/* Functions added for use with EIGBICG */
/* ------------------------------------ */
void BLAS_CGEEV(char *JOBVL, char *JOBVR, int *N, Complex_C *A, int *LDA, Complex_C *W, Complex_C *VL, int *LDVL,
                 Complex_C *VR, int *LDVR, Complex_C *WORK, int *LWORK, float *RWORK, int *INFO); 

void BLAS_ZGEEV(char *JOBVL, char *JOBVR, int *N, Complex_Z *A, int *LDA, Complex_Z *W, Complex_Z *VL, int *LDVL,
                 Complex_Z *VR, int *LDVR, Complex_Z *WORK, int *LWORK, double *RWORK, int *INFO); 

void BLAS_CGESV(int *N, int *NRHS, Complex_C *A, int *LDA, int *IPIV, Complex_C *B, int *LDB, int* INFO);

void BLAS_ZGESV(int *N, int *NRHS, Complex_Z *A, int *LDA, int *IPIV, Complex_Z *B, int *LDB, int* INFO);

void BLAS_CGETRF(int *M, int *N, Complex_C *A, int *LDA, int *IPIV, int *INFO);

void BLAS_CGETRS(char *TRANS, int *N, int *NRHS, Complex_C *A, int *LDA, int *IPIV, Complex_C *B, int *LDB, int *INFO);

void BLAS_ZGETRF(int *M, int *N, Complex_Z *A, int *LDA, int *IPIV, int *INFO);

void BLAS_ZGETRS(char *TRANS, int *N, int *NRHS, Complex_Z *A, int *LDA, int *IPIV, Complex_Z *B, int *LDB, int *INFO);

/*-------------------------------------------------------------------------------------------------------------------------*/








#ifdef NUM_ESSL
/* Use these in ESSL */
int dspev(int iopt, double *ap, double *w, double *z, int ldz, int n, double *aux, int naux);
int zhpev(int iopt, void *ap, double *w, void *z, int ldz, int n, void *aux, int naux);
#endif

/* END OF GENERAL NON CRAY ONES */
#else
/* These should be used for CRAY */

void BLAS_DCOPY(int *n, double *x, int *incx, double *y, int *incy);
void BLAS_DSWAP(int *n, double *x, int *incx, double *y, int *incy);
void BLAS_DGEMM(_fcd transa_fcd, _fcd transb_fcd, int *m, int *n, int *k, 
   double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, 
   double *c, int *ldc);
void BLAS_DSYMM(_fcd side_fcd, _fcd uplo_fcd, int *m, int *n, double *alpha, 
   double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void BLAS_DAXPY(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void BLAS_DGEMV(_fcd transa_fcd, int *m, int *n, double *alpha, double *a, int *lda, 
   double *x, int *incx, double *beta, double *y, int *incy);
double BLAS_DDOT(int *n, double *x, int *incx, double *y, int *incy);
double BLAS_DLAMCH(_fcd cmach_fcd);
void BLAS_DSCAL(int *n, double *alpha, double *x, int *incx);
void BLAS_DLARNV(int *idist, int *iseed, int *n, double *x);
void BLAS_DSYEV(_fcd jobz_fcd, _fcd uplo_fcd, int *n, double *a, int *lda, double *w,
   double *work, int *ldwork, int *info);

void BLAS_DSYTRF(_fcd uplo, int *n, double *a, int *lda, int *ipivot, double *work,
   int *ldwork, int *info);
void BLAS_DSYTRS(_fcd uplo, int *n, int *nrhs, double *a, int *lda, int *ipivot,
   double *b, int *ldb, int *info);

void   BLAS_ZCOPY(int *n, void *x, int *incx, void *y, int *incy);
void   BLAS_ZSWAP(int *n, void *x, int *incx, void *y, int *incy);
void   BLAS_ZGEMM(_fcd transa, _fcd transb, int *m, int *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   BLAS_ZHEMM(_fcd side, _fcd uplo, int *m, int *n, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   BLAS_ZAXPY(int *n, void *alpha, void *x, int *incx, void *y, int *incy);
void   BLAS_ZGEMV(_fcd transa, int *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
void   BLAS_ZSCAL(int *n, void *alpha, void *x, int *incx);
void   BLAS_ZLARNV(int *idist, int *iseed, int *n, void *x);
void   BLAS_ZHEEV(_fcd jobz, _fcd uplo, int *n, void *a, int *lda, double *w, void *work, int *ldwork, double *rwork, int *info);
void   ZDOTCSUB(void *dot, int *n, void *x, int *incx, void *y, int *incy);

void   BLAS_CPOTRF(_fcd UPLO, int *N, void *A, int *LDA, int *info );
void   BLAS_CPOTRS(_fcd UPLO, int *N, int *NRHS, Complex_C *A, int *LDA, B, int *LDB,int *INFO);
void   BLAS_ZHETRF(_fcd uplo, int *n, void *a, int *lda, int *ipivot, void *work, int *ldwork, int *info);
void   BLAS_ZHETRS(_fcd uplo, int *n, int *nrhs, void *a, int *lda, int *ipivot, void *b, int *ldb, int *info);

#endif /* else (if cray)*/

#ifdef Cplusplus
}
#endif /* Cplusplus */

#endif /* NUMERICAL_H */
