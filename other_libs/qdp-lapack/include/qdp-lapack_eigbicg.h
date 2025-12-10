// $Id: qdp-lapack_eigbicg.h,v 1.2 2009-10-21 21:03:39 kostas Exp $
/*******************************************************************************
    Eig-BiCG algorithm for solving ono-symmetric linear system Ax=b and 
    approximating few small eignvalues and eignvectors.
    
    Authors: Abdou M. Abdel-Rehim, Andreas Stathopolous , and Kostas Orginos
             andreas@cs.wm.edu

    Last updated: 10/14/2009
*******************************************************************************/ 




#ifndef __EIGBICG_C_H
#define __EIGBICG_C_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "qdp-lapack_restart_X.h"           
#include "qdp-lapack_biortho.h" 

#ifdef __cplusplus
extern "C" {
#endif 


/*******************************************************************************/
void Ceigbicg(int n, int lde, Complex_C *x, Complex_C *b, float *normb, 
              float tol, int maxit, char SRT_OPT, float epsi, int ConvTestOpt, 
              int *iter, float *reshist, int *flag, int plvl, Complex_C *work, 
              void (*matvec) (void *, void *, void *), 
              void (*mathvec)(void *, void *, void *), void *params, 
              float *AnormEst, int nev, Complex_C *evals, float *rnorms, 
              int v_max, Complex_C *VR, int LDVR, Complex_C *VL, int LDVL, 
              int esize, Complex_C *ework);

/*  Ceigbicg- Solve Ax=b using BICG for a non-Hermitian matrix A.
  	     At the same time compute few eigenvalues and vectors of A
             with smallest absolute value or smallest real part.
            
             Matrix A is non-Hermitian. It is accessed by 
	     matvec(x,y,params) and mathvec(x,y,params) where y,x are single 
             precision complex vectors. params is a void parameter where in 
             the user can encapsulate all needed information for 
             matvec and mathvec ops.

             In this version we store vectors in single precision but the global 
             sums are done in double precision and the result stored in a single 
             precision variable. The projection matrix and the bicg paramters are 
             treated as double precision.
           

  Parameters:
  ----------------------------------
  n        active problem size, The active part of A is an N-by-N matrix
  lde	   physical problem size. A is stored as (lde x n) array.
  ----------------------------------
            BICG-related parameters
  ----------------------------------
  x         (IN) the initial guess
            (OUT) the computed approximate solution
  b         (IN) the right hand side of the system
  tol       (IN) error tolerance r < tol*||b||
  maxit     (IN) maximum number of iterations
  SRT_OPT   (IN) 'M' means eigenvalues will be with smallest absolute value,
                 'R' means eigenvalues will be with smallest real part.

  epsi      (IN) thersold value used to determine if two numbers are conjugate
                 pairs within epsi. (see conjugate.h).

  ConvTestOpt (IN): option for how to determine convergence of the linear system
                    1  means norm(residual) < tol*norm(b)
                    2  means norm(residual) < MAX( tol*norm(b), MACHEPS*(AnormEst*norm(x)+norm(b)))
                    with MACHEPS=1e-7 for single precision version and 1e-16 for double precision version.

  reshist   (OUT) convergence history of residual norms for the iter iterations.

  iter      (OUT) number BICG iterations performed.
  flag      (OUT) exit status (see below)
  plvl      print info level (0: no output, 1/2/3 levels see below)
  work      (IN/OUT) work array. Must be of size 6*lde >= 6*n
  matvec    function that performs multiplication with matrix A
  mathvec   function that performs multiplication with matrix A^H
  params    (IN) parameters needed to do matrix-std::vector multiplications.
  AnormEst  (IN/OUT) input estimate of norm(A). If the absolute value 
                     largest eigenvalue of the v_max eigenvalues computed 
                     while restarting is larger, it will be replaced on output. 
  ----------------------------------
            eigen-related parameters
  ----------------------------------
  nev      (IN) number of eigenvalues to find
  evals    (OUT)nev Ritz values with smallest real part or smallest absolute value.
  v_max    (IN) maximum number of basis vectors used to approximate eigenvalues. 
  rnorms   (OUT)residual norms of the nev returned eigenpairs computed using:
                 rnorms[i] = norm(A*VR_i - evals_i*VR_i)/norm(VR_i).
  VR	   (OUT) the first (LDVR \times nev) contain the right Ritz vectors, std::vector by 
   	         std::vector. Users may then copy them to the desired data structure
  LDVR     (IN)  Leading dimesnion of right Ritz vectors.
  VL	   (OUT) the first (LDVL \times nev) contain the left Ritz vectors, std::vector by 
   	         std::vector. Users may then copy them to the desired data structure.
  LDVL     (IN) Leading dimension of left Ritz vectors.

  esize	   (IN) size of ework, the eigenwork space: 
  		 2*lde + 2*nev =< esize =< 2(nev+1)*lde 
  ework    (IN/OUT)temp work space of size esize
  ----------------------------------

  On exit, if flag is

   0 then BICG converged to the desired tolerance TOL within MAXIT iterations 
   	convergence test:  norm(residual) < tol*norm(b)

   1 then BICG iterated MAXIT times but did not converge.

   2 then one of the scalar quantities computed during BICG was zero

   3 then error occured while using LAPCK to solve the projected small eigenvalue problem
   ----------------------------------
   PLVL    0 no output
   	 >=1 prints BICG linear system info on exit (num its/ flag)
	 >=2 prints linear system residuals at every iteration 
	 >=3 prints final eigenvalues and their residual norm as well as
             the norm of the left eigenstd::vector and the cos of the angle between left and
             right eigenvectors "cangle=vl'*vr/norm(vl)/norm(vr)". 
	 >=4 prints eigenvalues at every restart and their residual norm as well as
             the norm of the left eigenstd::vector and the cos of the angle between left and
             right eigenvectors "cangle=vl'*vr/norm(vl)/norm(vr)". 

*****************************************************************************************************/

void Ceigbicg_ver2(int n, int lde, Complex_C *x, Complex_C *b, float *normb, 
              float tol, int maxit, char SRT_OPT, float epsi, int ConvTestOpt, 
              int *iter, float *reshist, int *flag, int plvl, Complex_C *work, 
              void (*matvec) (void *, void *, void *), 
              void (*mathvec)(void *, void *, void *), void *params, 
              void (*dmatvec)(void *, void *, void *), void *dparams, 
              float *AnormEst, int nev, Complex_C *evals, float *rnorms, 
              int v_max, Complex_Z *VR, int LDVR, Complex_Z *VL, int LDVL, 
              int esize, Complex_Z *ework);




/*****************************************************************************************************/
void Zeigbicg(int n, int lde, Complex_Z *x, Complex_Z *b, double *normb, 
              double tol, int maxit, char SRT_OPT, double epsi, int ConvTestOpt, int *iter, 
              double *reshist, int *flag, int plvl, Complex_Z *work, 
              void (*matvec) (void *, void *, void *), 
              void (*mathvec)(void *, void *, void *), void *params, 
              double *AnormEst, int nev, Complex_Z *evals, double *rnorms, 
              int v_max, Complex_Z *VR, int LDVR, Complex_Z *VL, int LDVL, 
              int esize, Complex_Z *ework);


/* In this version, all variables are treated in double precision*/ 

/****************************************************************************************************/
static void CdisplayInfo(float tol, int maxit, int flag, int iter, float resnorm);
/* displayInfo- prints information about the linear system solution with Ceigbicg on exit
                to the stdout.
   tol     (IN): The requested tolerance of convergence.
   maxit   (IN): Maximum number of iterations allowed.
   flag    (IN): The exit status flag of eigbicg.
   iter    (IN): Number of iterations performed by eigbicg on exit.
   resnorm (IN): The residual norm |b-Ax| on exit. */
/***************************************************************************************************/
static void ZdisplayInfo(double tol, int maxit, int flag, int iter, double resnorm);

// double version for displayInfo



/*---------------------------------------------------------------------------------------------------------------*/
void CcomputeResNorm( Complex_C *xr, Complex_C *xl, Complex_C *lambda, float *rnorm, int n, Complex_C *Res, 
                      float *xlnorm, float *xrnorm, Complex_C *cangle, 
                      void (*matvec)(void *, void *, void *), void *params);



/* computeResNorm- Given a right eigenstd::vector xr, it computes its Ritz value lambda=xr'*A*xr/norm(xr), its
   residual Res=Axr-lambda*xr, and its residual norm rnorm.
   Requires matvec and the parameters it uses params .

   xr (IN): right eigenstd::vector
   xl (IN): corresponding left eigenstd::vector.
   lambda (OUT): Ritz value xr'*A*xr /(xr'*xr).
   rnorm  (OUT): norm(A*xr-lambda*xr).
   Res    (OUT): eigenvalue residual A*xr-lambda*xr.
   xlnorm (OUT): norm(xl)
   xrnorm (OUT): norm(xr)
   cangle (OUT): angle between left and right eigenvectors cangle=xl'*xr/norm(xl)/norm(xr).
   matvec (IN) : function that performs matrix std::vector product.
   params (IN) : parameters needed by matvec.
-------------------------------------------------------------------------------------------------------------------*/
void ZcomputeResNorm( Complex_Z *xr, Complex_Z *xl, Complex_Z *lambda, double *rnorm, int n, Complex_Z *Res, 
                      double *xlnorm, double *xrnorm, Complex_Z *cangle, 
                      void (*matvec)(void *, void *, void *), void *params);


// double version of computeResNorm



#ifdef __cplusplus
};
#endif 
#endif
