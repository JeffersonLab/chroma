// $Id: qdp-lapack_IncrEigbicg.h,v 1.1 2009-10-21 20:50:56 kostas Exp $
/**************************************************************************************
  Incremental eigbicg for solving non-symmetric systems with multiple right-hand sides
  Authors: Abdou Abdel-Rehim, Andreas Stathopolous, and Kostas Orginos
           andreas@cs.wm.edu
  Last updated : 10/19/2009
***************************************************************************************/

#ifndef  __INCREIGBICG_C_H
#define  __INCREIGBICG_C_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "qdp-lapack_eigbicg.h"
#include "qdp-lapack_wtime.h"
#include "bicgstab.h"

#ifdef __cplusplus
extern "C" {
#endif 

/***********   SINGLE PRECISION VERSION WITH DOUBLE PRECISION GLOBAL SUM   **************/

void IncrEigbicg_C(
        int n, int lde,        /* (IN) n dim of matrix A, lde leading dim of A and vecs */
	int nrhs, 		/* (IN) The number of right hand sides to solve */
	Complex_C *X,   	/* (IN/OUT) The nrhs solutions */
	Complex_C *B,   	/* (IN) The nrhs right hand sides */
	int *ncurEvals,         /* (IN/OUT) Num of eigenpairs already in evecsl,evecsr,evals */
	int ldh,  		/* (IN) Max eigenpairs that can be stored in evecs */
                    
                                /* evecsl(ldvl*ldh), evecsr(ldvr*ldh), and H(ldh*ldh) has to be allocated by the calling program*/
	Complex_C *evecsl, 	/* (IN/OUT) The left eigenstd::vector basis augmented by eigbicg with ldvl leading dimension */
        Complex_C *evecsr, 	/* (IN/OUT) The right eigenstd::vector basis augmented by eigbicg with ldvr leading dimension */
	Complex_C *evals,       /* (OUT) The eigenvalues augmented by eigbicg */
	Complex_C *H,		/* (IN/OUT) The ncurEvals^2 matrix: H=evecsl'*A*evecsr; 
                                   H has ldh^2 storage space with ncurEvals^2 active part with ldh leading dimension*/
	                        /* -------------------------------------------- */
        void (*matvec) (void *, void *, void *),  /* (IN) Matvec operator */
	void (*mathvec)(void *, void *, void *),  /* (IN) Mat^H vec operator */
	void *params,  	        /* (IN) parameters to be passed to operators */
        float *AnormEst,       /*(IN/OUT) estimate of the norm of A. */
				/* ---------------------------------------- */
				
			        /* arrays work,VL,VR and ework will be allocated inside IncrEigbicg*/
	Complex_C *work,        /* work array for BICGCG etc of dimension  6*lde */
	Complex_C *VL, 		/*  (IN)work array for left  eigenstd::vector search basis of dimension ldvl*v_max */
        int       ldvl,         /* (IN) Leading dimension of left vectors */
	Complex_C *VR, 		/* (IN) work array for right eigenstd::vector search basis of dimension ldvr*v_max */
        int       ldvr,         /* (IN) Leading dimension of right vectors*/
	Complex_C *ework, 	/* (IN) work array. If esize < 2*lde+2*nev allocate with appropriate esize*/
	int esize, 		/* (IN) 2*lde+2*nev <= esize <= 2*(nev+1)*lde */
				/* ---------------------------------------- */
	float tol, 		/* (IN) BICG and BICGSTAB linear system tolerance |res| < tol*b */
	float *restartTol, 	/* (IN) Restart BICGSTAB when |res|<restartTol*|b-Ax(0)| */
	int maxit, 		/* (IN) Maximum number of BICG ind BICGSTAB terations per system */
        char SRT_OPT,          /* (IN) option for sorting eigenvalues computed by eigbicg
                                   'M' for smallest magnitude and 'R' for smallest real part*/
        float epsi,            /* (IN) threshold used to check if two eignvalues are conjugate pairs:
                                   if( imag(x)*imag(y) < 0  and
                                   abs(imag(x)+imag(y))/abs(x+y) < epsi and
                                   abs(real(x)-real(y)) / abs(x+y) < epsi )
                                   then x and y are considered to be conjugate, otherwise they are not.*/
        int ConvTestOpt,       /* (IN)option for how to determine convergence of the linear system
                                       1  means norm(residual) < tol*norm(b)
                                       2  means norm(residual) < MAX( tol*norm(b), MACHEPS*(AnormEst*norm(x)+norm(b)))
                                       with MACHEPS=1e-7 for single precision version and 1e-16 for double precision version.*/
	int plvl, 		/* (IN) Printing level [0..4] */
	int nev, 		/* (IN) Number of eigenvalues to target in each BICG */
	int v_max,		/* (IN) Maximum number of n-dim vectors in VL,VR */
	FILE *outputFile);      /* (IN) File to write reports. Could be STDERR */

/***************************************************************************************************/

/* init_BICG_C computes and improved initial guess xinit by performing a left-right projection over a 
   given spaces vl and vr. 
   Let xinit be the initial guess for the system A*x=b and that we have nvecs left and right 
   biorthogonal spaces vl and vr. An improved initial guess is given by xinit_new = xinit + vr *sol such
   that sol is given by the condition that the new residual is orthogonal to vl. 
   res_new = b - A*xinit_new and vl'*res_new=0. This gives:
   vl'*b - vl'*A*(xinit+vr*sol) = 0, vl'*(b-A*xinit) - vl'*A*vr*sol =0 , vl'*res - H*sol =0.
   So, sol is given by the solution sol = inv(H)*vl'*res and
   xinit_new= xinit + vr*(inv(H)*vl'*res) */
void init_BICG_C(
               Complex_C *vl,    /*(IN)left nvecs vectors */ 
               int ldvl,         /*(IN)leading dimension of vl, n is the active part*/
               Complex_C *vr,    /*(IN)right nvecs vectors such that vl'*vr=I */
               int ldvr,         /*(IN)leading dimension of vr, n is the active part*/
               int nvecs,        /*(IN)number of left and right vectors used for the projection*/
               Complex_C *xinit, /*(IN/OUT) on input is the current initial guess, 
                                            on output is the projected intial guess. 
                                            Leading dimension lde and active dimension n*/
               Complex_C *b,     /*(IN) right-hand side with leading dimension lde and active dimension n*/
               int lde,          /*(IN) leading dimension of A, xinit and b */
               int n,            /*(IN) active A dimension*/
               Complex_C *H,     /*(IN) projection matrix vl'*A*vr of dimension (ldh*nvecs) with nvecs active part*/
               int ldH,          /*(IN) leading dimension of H */
               int *IPIV,        /*(IN/OUT) nvecs dimension integer array used by the lapack linear solver */
               Complex_C *work,  /*(IN/OUT) work array of dimension lde+nvecs */
               void (*matvec) (void *, void *, void *), /* (IN) matrix std::vector operator */
               void *params,     /*(IN) paramters needed by matvec */
               int *info         /*(OUT) output status of LAPACK CGESV linear solver :
                                         0 means sucessful exit
                                         >0 if info=i then U(i,i) (LU factor) was zero.
                                         <0, if info=-i, it means the i-th argument had an illegal value*/ 
               );

/*********************************************************************************************************************/

/* LRD_BICGSTAB_C is left-right deflated restarted bicgstab. Given "k" left and right biorthogonal subspaces "vl" and "vr"
   , the projection matrix "H=vl^H*A*vr", restarting tolerance "DefTol" and convergence tolerance "tol", this function
   will perform deflation using left-right projection with init_BICG and then use bicgstab. When the system reaches
   DefTol convergence level, the projection step is performed one more time and then bicgstab is called again. */

void LRD_BICGSTAB_C(
               Complex_C *vl,    /*(IN)left nvecs vectors */ 
               int ldvl,         /*(IN)leading dimension of vl, n is the active part*/
               Complex_C *vr,    /*(IN)right nvecs vectors such that vl'*vr=I */
               int ldvr,         /*(IN)leading dimension of vr, n is the active part*/
               int nvecs,        /*(IN)number of left and right vectors used for the projection*/
               Complex_C *x    , /*(IN/OUT) on input is the initial guess, 
                                            on output is the solution. 
                                            Leading dimension lde and active dimension n*/
               Complex_C *b,     /*(IN) right-hand side with leading dimension lde and active dimension n*/
               int lde,          /*(IN) leading dimension of A, xinit and b */
               int n,            /*(IN) active A dimension*/
               Complex_C *H,     /*(IN) projection matrix vl'*A*vr of dimension (ldh*nvecs) with nvecs active part*/
               int ldH,          /*(IN) leading dimension of H */
               int *IPIV,        /*(IN/OUT) permutations used in LU factorization */
               Complex_C *work,  /*(IN/OUT) work array of dimension 6*lde*/
               void (*matvec) (void *, void *, void *), /* (IN) matrix std::vector operator */
               void *params,     /*(IN) paramters needed by matvec */
               float AnormEst,   /* (IN) estimate value of the norm of A */
               int maxiter,      /*maximum number of iterations that could be used in bicgstab*/
               float DefTol,     /*restarting tolerance*/
               float tol,        /* convergence tolerance of the linear system. At convergence norm(b-A*x) < tol*norm(b) */
               int ConvTestOpt,       /* (IN)option for how to determine convergence of the linear system
                                       1  means norm(residual) < tol*norm(b)
                                       2  means norm(residual) < MAX( tol*norm(b), MACHEPS*(AnormEst*norm(x)+norm(b)))
                                       with MACHEPS=1e-7 for single precision version and 1e-16 for double precision version.*/  
               FILE *outputFile, /* output file stream, could be stout or stderr*/
               int *info         /*(OUT) output status */ 
               );

/***************************************************************************************************************************/

/* ComputeFinalEvecs_C will compute evals and evecs from a left and right basis evecsl,evecsr which is biorthogonal and 
                     the projection matrix H. It computes the Ritz values and residual norms.*/

void ComputeFinalEvecs_C
    ( int nvecs, int n, Complex_C *evecsl, int ldvl, Complex_C *evecsr, int ldvr, Complex_C *H, int ldh, char SRT_OPT, 
      float epsi, Complex_C *evals, float *ernorms, float *xlnorms, float *xrnorms, Complex_C *angles, 
      void (*matvec)(void *, void *, void *), void *params, Complex_C *work, int worksize);

/* given nvecs left and right basis vectors evecsl and evecsr as well as the projection matrix
   H = evecsl' *A * evecsr, this function computes Ritz values evals, Ritz vectors evecsl, and evecsr
   and the eigenvalues resudal norms ernorms. It also returns the norms and angles of the left and right
   vectors in xlnorms, xrnorms, and angles.  

   
   nvecs(IN) : number of bais vectors.
   evecsl(IN/OUT): (IN) nvecs left basis, (OUT) nvecs left eigenvectors
   ldvl(IN)  : leading dimensions of left vectors.
   evcsr(IN/OUT): (IN) nvecs right basis, (OUT) nvecs right eigenvectors
   ldvr(IN): leading dimension of right vectors.
   H (IN)  : projection matrix. UN-destroyed on output.
   ldh(IN) : leading dimension of H.
   SRT_OPT (IN): sorting criteria of the eigenvectors
                 'M' means sort in ascending order w.r.t. absolute value
                 'R' means sort in ascending order w.r.t real part.
   float epsi,   (IN) threshold used to check if two eignvalues are conjugate pairs:
                 if( imag(x)*imag(y) < 0  and
                 abs(imag(x)+imag(y))/abs(x+y) < epsi and
                 abs(real(x)-real(y)) / abs(x+y) < epsi )
                 then x and y are considered to be conjugate, otherwise they are not.
   evals(OUT): nvecs Ritz values ordered according to SRT_OPT.
   ernorms(OUT): norms of eigenvalues residuals=norms(A*v-lambda*v)/norm(v)
   xlnorms(OUT): norms of left eigenvectors.
   xrnorms(OUT): norms of right eigenvectors.
   angles(OUT) : vl'*vr / norm(vl)/norm(vr).
   matvec (IN) : function that performs matrix std::vector product.
   params (IN) : parameters needed by matvec.
-------------------------------------------------------------------------------------------------------------------*/


/******************************DOUBLE PRECISION VERSION**************************************************/


void IncrEigbicg_Z(
        int n, int lde,        /* (IN) n dim of matrix A, lde leading dim of A and vecs */
	int nrhs, 		/* (IN) The number of right hand sides to solve */
	Complex_Z *X,   	/* (IN/OUT) The nrhs solutions */
	Complex_Z *B,   	/* (IN) The nrhs right hand sides */
	int *ncurEvals,         /* (IN/OUT) Num of eigenpairs already in evecsl,evecsr,evals */
	int ldh,  		/* (IN) Max eigenpairs that can be stored in evecs */
                    
                                /* evecsl(ldvl*ldh), evecsr(ldvr*ldh), and H(ldh*ldh) has to be allocated by the calling program*/
	Complex_Z *evecsl, 	/* (IN/OUT) The left eigenstd::vector basis augmented by eigbicg with ldvl leading dimension */
        Complex_Z *evecsr, 	/* (IN/OUT) The right eigenstd::vector basis augmented by eigbicg with ldvr leading dimension */
	Complex_Z *evals,       /* (OUT) The eigenvalues augmented by eigbicg */
	Complex_Z *H,		/* (IN/OUT) The ncurEvals^2 matrix: H=evecsl'*A*evecsr; 
                                   H has ldh^2 storage space with ncurEvals^2 active part with ldh leading dimension*/
	                        /* -------------------------------------------- */
        void (*matvec) (void *, void *, void *),  /* (IN) Matvec operator */
	void (*mathvec)(void *, void *, void *),  /* (IN) Mat^H vec operator */
	void *params,  	        /* (IN) parameters to be passed to operators */
        double *AnormEst,       /*(IN/OUT) estimate of the norm of A. */
				/* ---------------------------------------- */
				
			        /* arrays work,VL,VR and ework will be allocated inside IncrEigbicg*/
	Complex_Z *work,        /* work array for BICGCG etc of dimension  6*lde */
	Complex_Z *VL, 		/*  (IN)work array for left  eigenstd::vector search basis of dimension ldvl*v_max */
        int       ldvl,         /* (IN) Leading dimension of left vectors */
	Complex_Z *VR, 		/* (IN) work array for right eigenstd::vector search basis of dimension ldvr*v_max */
        int       ldvr,         /* (IN) Leading dimension of right vectors*/
	Complex_Z *ework, 	/* (IN) work array. If esize < 2*lde+2*nev allocate with appropriate esize*/
	int esize, 		/* (IN) 2*lde+2*nev <= esize <= 2*(nev+1)*lde */
				/* ---------------------------------------- */
	double tol, 		/* (IN) BICG and BICGSTAB linear system tolerance |res| < tol*b */
	double *restartTol, 	/* (IN) Restart BICGSTAB when |res|<restartTol*|b-Ax(0)| */
	int maxit, 		/* (IN) Maximum number of BICG ind BICGSTAB terations per system */
        char SRT_OPT,          /* (IN) option for sorting eigenvalues computed by eigbicg
                                   'M' for smallest magnitude and 'R' for smallest real part*/
        double epsi,            /* (IN) threshold used to check if two eignvalues are conjugate pairs:
                                   if( imag(x)*imag(y) < 0  and
                                   abs(imag(x)+imag(y))/abs(x+y) < epsi and
                                   abs(real(x)-real(y)) / abs(x+y) < epsi )
                                   then x and y are considered to be conjugate, otherwise they are not.*/
        int ConvTestOpt,       /* (IN)option for how to determine convergence of the linear system
                                       1  means norm(residual) < tol*norm(b)
                                       2  means norm(residual) < MAX( tol*norm(b), MACHEPS*(AnormEst*norm(x)+norm(b)))
                                       with MACHEPS=1e-7 for single precision version and 1e-16 for double precision version.*/
	int plvl, 		/* (IN) Printing level [0..4] */
	int nev, 		/* (IN) Number of eigenvalues to target in each BICG */
	int v_max,		/* (IN) Maximum number of n-dim vectors in VL,VR */
	FILE *outputFile);      /* (IN) File to write reports. Could be STDERR */


/**************************************************************************************************************/

/* init_BICG_Z computes and improved initial guess xinit by performing a left-right projection over a 
   given spaces vl and vr. 
   Let xinit be the initial guess for the system A*x=b and that we have nvecs left and right 
   biorthogonal spaces vl and vr. An improved initial guess is given by xinit_new = xinit + vr *sol such
   that sol is given by the condition that the new residual is orthogonal to vl. 
   res_new = b - A*xinit_new and vl'*res_new=0. This gives:
   vl'*b - vl'*A*(xinit+vr*sol) = 0, vl'*(b-A*xinit) - vl'*A*vr*sol =0 , vl'*res - H*sol =0.
   So, sol is given by the solution sol = inv(H)*vl'*res and
   xinit_new= xinit + vr*(inv(H)*vl'*res) */
void init_BICG_Z(
               Complex_Z *vl,    /*(IN)left nvecs vectors */ 
               int ldvl,         /*(IN)leading dimension of vl, n is the active part*/
               Complex_Z *vr,    /*(IN)right nvecs vectors such that vl'*vr=I */
               int ldvr,         /*(IN)leading dimension of vr, n is the active part*/
               int nvecs,        /*(IN)number of left and right vectors used for the projection*/
               Complex_Z *xinit, /*(IN/OUT) on input is the current initial guess, 
                                            on output is the projected intial guess. 
                                            Leading dimension lde and active dimension n*/
               Complex_Z *b,     /*(IN) right-hand side with leading dimension lde and active dimension n*/
               int lde,          /*(IN) leading dimension of A, xinit and b */
               int n,            /*(IN) active A dimension*/
               Complex_Z *H,     /*(IN) projection matrix vl'*A*vr of dimension (ldh*nvecs) with nvecs active part*/
               int ldH,          /*(IN) leading dimension of H */
               int *IPIV,        /*(IN/OUT) nvecs dimension integer array used by the lapack linear solver */
               Complex_Z *work,  /*(IN/OUT) work array of dimension lde+nvecs */
               void (*matvec) (void *, void *, void *), /* (IN) matrix std::vector operator */
               void *params,     /*(IN) paramters needed by matvec */
               int *info         /*(OUT) output status of LAPACK CGESV linear solver :
                                         0 means sucessful exit
                                         >0 if info=i then U(i,i) (LU factor) was zero.
                                         <0, if info=-i, it means the i-th argument had an illegal value*/ 
               );

/*********************************************************************************************************************/

/* LRD_BICGSTAB_Z is left-right deflated restarted bicgstab. Given "k" left and right biorthogonal subspaces "vl" and "vr"
   , the projection matrix "H=vl^H*A*vr", restarting tolerance "DefTol" and convergence tolerance "tol", this function
   will perform deflation using left-right projection with init_BICG and then use bicgstab. When the system reaches
   DefTol convergence level, the projection step is performed one more time and then bicgstab is called again. */

void LRD_BICGSTAB_Z(
               Complex_Z *vl,    /*(IN)left nvecs vectors */ 
               int ldvl,         /*(IN)leading dimension of vl, n is the active part*/
               Complex_Z *vr,    /*(IN)right nvecs vectors such that vl'*vr=I */
               int ldvr,         /*(IN)leading dimension of vr, n is the active part*/
               int nvecs,        /*(IN)number of left and right vectors used for the projection*/
               Complex_Z *x    , /*(IN/OUT) on input is the initial guess, 
                                            on output is the solution. 
                                            Leading dimension lde and active dimension n*/
               Complex_Z *b,     /*(IN) right-hand side with leading dimension lde and active dimension n*/
               int lde,          /*(IN) leading dimension of A, xinit and b */
               int n,            /*(IN) active A dimension*/
               Complex_Z *H,     /*(IN) projection matrix vl'*A*vr of dimension (ldh*nvecs) with nvecs active part*/
               int ldH,          /*(IN) leading dimension of H */
               int *IPIV,        /*(IN/OUT) permutations used in LU factorization */
               Complex_Z *work,  /*(IN/OUT) work array of dimension 6*lde*/
               void (*matvec) (void *, void *, void *), /* (IN) matrix std::vector operator */
               void *params,     /*(IN) paramters needed by matvec */
               double AnormEst,   /* (IN) estimate value of the norm of A */
               int maxiter,      /*maximum number of iterations that could be used in bicgstab*/
               double DefTol,     /*restarting tolerance*/
               double tol,        /* convergence tolerance of the linear system. At convergence norm(b-A*x) < tol*norm(b) */
               int ConvTestOpt,       /* (IN)option for how to determine convergence of the linear system
                                       1  means norm(residual) < tol*norm(b)
                                       2  means norm(residual) < MAX( tol*norm(b), MACHEPS*(AnormEst*norm(x)+norm(b)))
                                       with MACHEPS=1e-7 for single precision version and 1e-16 for double precision version.*/  
               FILE *outputFile, /* output file stream, could be stout or stderr*/
               int *info         /*(OUT) output status */ 
               );

/***************************************************************************************************************************/

/* ComputeFinalEvecs_Z will compute evals and evecs from a left and right basis evecsl,evecsr which is biorthogonal and 
                     the projection matrix H. It computes the Ritz values and residual norms.*/

void ComputeFinalEvecs_Z
    ( int nvecs, int n, Complex_Z *evecsl, int ldvl, Complex_Z *evecsr, int ldvr, Complex_Z *H, int ldh, char SRT_OPT, 
      double epsi, Complex_Z *evals, double *ernorms, double *xlnorms, double *xrnorms, Complex_Z *angles, 
      void (*matvec)(void *, void *, void *), void *params, Complex_Z *work, int worksize);

/* given nvecs left and right basis vectors evecsl and evecsr as well as the projection matrix
   H = evecsl' *A * evecsr, this function computes Ritz values evals, Ritz vectors evecsl, and evecsr
   and the eigenvalues resudal norms ernorms. It also returns the norms and angles of the left and right
   vectors in xlnorms, xrnorms, and angles.  

   
   nvecs(IN) : number of bais vectors.
   evecsl(IN/OUT): (IN) nvecs left basis, (OUT) nvecs left eigenvectors
   ldvl(IN)  : leading dimensions of left vectors.
   evcsr(IN/OUT): (IN) nvecs right basis, (OUT) nvecs right eigenvectors
   ldvr(IN): leading dimension of right vectors.
   H (IN)  : projection matrix. UN-destroyed on output.
   ldh(IN) : leading dimension of H.
   SRT_OPT (IN): sorting criteria of the eigenvectors
                 'M' means sort in ascending order w.r.t. absolute value
                 'R' means sort in ascending order w.r.t real part.
   float epsi,   (IN) threshold used to check if two eignvalues are conjugate pairs:
                 if( imag(x)*imag(y) < 0  and
                 abs(imag(x)+imag(y))/abs(x+y) < epsi and
                 abs(real(x)-real(y)) / abs(x+y) < epsi )
                 then x and y are considered to be conjugate, otherwise they are not.
   evals(OUT): nvecs Ritz values ordered according to SRT_OPT.
   ernorms(OUT): norms of eigenvalues residuals=norms(A*v-lambda*v)/norm(v)
   xlnorms(OUT): norms of left eigenvectors.
   xrnorms(OUT): norms of right eigenvectors.
   angles(OUT) : vl'*vr / norm(vl)/norm(vr).
   matvec (IN) : function that performs matrix std::vector product.
   params (IN) : parameters needed by matvec.
-------------------------------------------------------------------------------------------------------------------*/
                 




#ifdef __cplusplus
};
#endif 
#endif
