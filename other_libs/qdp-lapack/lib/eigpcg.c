/* $Id: eigpcg.c,v 1.9 2009-10-21 20:50:57 kostas Exp $ */
/*-------------------------------------------------------------------------

  EIGPCG   Solve Ax=b by the preconditioned conjugate gradient method.
  	   At the same time compute smallest abs eigenvalues/vectors of A.
           Refs: Golub and Van Loan
	         Stathopoulos and Orginos
           
           Matrix A is Hermitian positive definite. It is accessed by 
	   matvec(x,y,params) where y,x are single precision complex vectors.
	   params is a void parameter where in the user can encapsulate 
	   all needed information for matvec and preconditioning ops.

  Parameters:
  ----------------------------------
  N        active problem size, The active part of A is an N-by-N matrix
  LDE	   physical problem size. A and vectors are stored in LDE size
  ----------------------------------
            CG-related parameters
  ----------------------------------
  X         (IN) the initial guess
            (OUT) the computed approximate solution
  B         (IN) the right hand side of the system
  NORMB     (IN/OUT) ||B|| is computed. On input, if flag==3, NORMB=||B||
  TOL       (IN) error tolerance r < tol*||b||
  RESTARTTOL(IN) restart CG when r < tol*||b-Ax0||
  MAXIT     (IN) maximum number of iterations
  RESHIST   (OUT) convergence history of residual norms (*of size MAXIT*)
  ITER      (IN/OUT) number PCG iterations performed in previous restarts (IN) 
  		     and previous+current iterations total (OUT)
  FLAG      (OUT) exit status (see below)
  PLVL      print info level (0: no output, 1/2/3 levels see below)
  WORK      (IN/OUT) work array. Must be of size 4*LDE >= 4*N
  MATVEC    function that performs multiplication with matrix A
  PRECON    function that applies preconditioner (can be NULL)
  ----------------------------------
            eigen-related parameters
  ----------------------------------
  NEV      (IN) number of eigenvalues to find
  M_MAX    (IN) maximum number of basis vectors
  EVALS    (OUT) nev Ritz values from smallest to largest
  RNORMS   (OUT) residual norms of the nev returned eigenpairs
  V	   (IN) the basis vectors  (N \times v_max)
           (OUT) the first (N\times nev) contain the Ritz vectors, std::vector by 
   	         std::vector. Users may then copy them to the desired data structure
  ESIZE	   (IN) size of ework, the eigenwork space: the more the better
  		N+2*nev <= esize <= (2*nev+1)*N
  EWORK    temp work space of size esize
  ----------------------------------

  On exit, if FLAG is

   0 then PCG converged to the desired tolerance TOL within MAXIT iterations 
   	convergence test:  norm(preconditioned residual) < TOL*norm(b)

   1 then PCG iterated MAXIT times but did not converge.

   2 then one of the scalar quantities computed during PCG was zero

   3 then PCG stopped because the restarting tolerance (related to initCG)
     is satisfied: i.e.,  norm(r(k)) < RestartTol*norm(r(1))
     Note: r(1)=b-Ax0 could be << r(0)=norm(b) 

   ----------------------------------
   PLVL    0 no output
   	 >=1 prints CG linear system info on exit (num its/ flag)
	 >=2 prints linear system residuals at every iteration 
	 >=3 prints norms of eigenvalue residuals at every restart using
	     the Lanczos estimates (cheap but inaccurate close to convergence)
	 >=4 prints norms of eigenvalue residuals by computing Ax-lx.
	     Very expensive -- only for debugging purposes. Also, assumes:
	     		ESIZE > 2*LDE + 2*nev
   ----------------------------------

   BLAS_<name> are macros that in numerical.h are defined to the
   appropriate way to call the Fortran BLAS <name> function. Typically 
   this involves an underscore at the end, eg. BLAS_saxpy -> saxpy_

   wrap_<name> is used for two functions. dot products and gemv for large
   parallel distributed vectors. The wrappers include a globalSum of the result

--------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // is this needed?
#include "qdp-lapack_eigpcg.h"

/* DEBUG FUNCTION **********************************************************/
void showSubMat(Complex_C *H, int lde, int rows, int cols){
	int i,j;
	printf(" Dim %dx%d\n",rows,cols);
	for (i=0;i<rows;i++) {
            for (j=0;j<cols;j++)
	       printf("(%g %g)",H[j*lde+i].r,H[j*lde+i].i);
	    printf("\n");
	}
}

/* EIGPCG *******************************************************************/
void eigpcg(int n, int lde, Complex_C *x, Complex_C *b, 
	 float *normb, float tol, float restartTol, int maxit,
	 int *iter, float *reshist, int *flag, int plvl, Complex_C *work,
	 void (*matvec)(void *, void *, void *),
	 void (*precon)(void *, void *, void *),
	 void *params,
	 int nev, float *evals, float *rnorms,
  	 int v_max, Complex_C *V, int esize, Complex_C *ework
	 )
{
  float tolb;			/* requested tolerance for residual */
  double alpha, beta;		/* CG scalars */
  double rho, rhoprev;
  double pAp;
  int it;			/* current iteration number */
  int i, j;			/* loop variables */
  int zs, cs, ds, tmpsize;

  Complex_C *r, *z, *p, *Ap;	/* ptrs in work for PCG vectors */

  Complex_Z tempz;		/* double precision complex temp var */
  Complex_C tempc;		/* single precision complex temp var */
  double    tempd;		/* double temp var */
  float     tempf;		/* float temp var */
  int tempi;			/* int temp var */
  int ONE = 1;			/* var for passing 1 into BLAS routines */
  /*----------------------------------------------------------------------
         Eigen variables and setup    
    ----------------------------------------------------------------------*/
  /* Some constants */
  Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
  char cR = 'R'; char cL = 'L'; char cN ='N'; 
  char cV = 'V'; char cU = 'U'; char cC ='C';
  double betaprev,
	 alphaprev;     /* remember the previous iterations scalars */
  int v_size;		/* tracks the size of V */
  int lwork = 2*v_max;  /* the size of zwork */
  Complex_C *Ap_prev;	/* ptr to first ework std::vector for storing Ap */
  Complex_C *tmp_work;  /* ptr to ework(N) for use in restarting */
  Complex_C *Coef;      /* single precision copy of restart coefs (v_max*2nev)*/
  Complex_Z *H;         /* the V'AV projection matrix */
  Complex_Z *Hevecs;    /* the eigenvectors of H */
  Complex_Z *Hevecsold; /* the eigenvectors of H(m-1,m-1) */
  double    *Hevals;    /* the eigenvalues of H */
  double    *Hevalsold; /* the eigenvalues of H(m-1,m-1) */
  Complex_Z *TAU;	/* The orthogonalization coefficients */
  Complex_Z *zwork;     /* double complex work array needed by zheev */
  double *rwork;        /* double work array needed by zheev */

  zs = sizeof(Complex_Z); 
  cs = sizeof(Complex_C); 
  ds = sizeof(double);

  if ((H = calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate H\n");
  if ((Hevecs = calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecs\n");
  if ((Hevecsold = calloc(v_max*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevecsold\n");
  if ((Hevals = calloc(v_max, ds)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevals\n");
  if ((Hevalsold = calloc(v_max-1, ds)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Hevalsold\n");
  if ((TAU = calloc(2*nev, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate TAU\n");
  if ((zwork = calloc(2*v_max, zs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate zwork\n");
  if ((rwork = calloc(3*v_max, ds)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate rwork\n");
  if ((Coef = calloc(2*nev*v_max, cs)) == NULL) 
     fprintf(stderr, "ERROR Could not allocate Coef\n");

  /* setup pointers into ework */
  Ap_prev  = ework;
  tmp_work = ework+n;
  tmpsize = esize - n; /* leftover in tmp_work */

  /*----------------------------------------------------------------------*/

  /* setup pointers into work */
  r = work;
  z = work + lde;
  p = work + 2*lde;
  Ap = work + 3*lde;

  if(precon==NULL) z=r ;

  /*--------------------------------------------------------------------
     Initialization phase 
  --------------------------------------------------------------------*/
  if (*flag != 3) {
     /* If flag == 3, the eigCG is called after restart with the same b 
      * whose norm is already known in normb, so no need for these    */
     tempc = wrap_cdot(&n, b, &ONE, b, &ONE, params); /* Norm of rhs, b */
     *normb = sqrt(tempc.r);

     /* If right hand side is zero return zero solution. ITER stays the same */
     if (*normb == 0.0) {		
       for (i=0; i<n; i++) { x[i].r=0.0; x[i].i=0.0;}
       *flag = 0;		
       reshist[0] = 0.0;
       if (plvl) displayInfo(tol,maxit,*flag,*iter,reshist[0]);
       return;
     }
  }
  
  /* Set up for the method */
  *flag = 1;
  tolb = tol * (*normb);	/* Relative to b tolerance */

  /* Zero-th residual: r = b - A*x  */
  matvec(x, r, params);			
  for (i = 0; i < n; i ++) {
      r[i].r = b[i].r - r[i].r;
      r[i].i = b[i].i - r[i].i;
  }
  
  rho = 0.0;
  alpha = 1.0;
  beta = 0.0;
  v_size = 0;

  /*--------------------------------------------------------------------
     main CG loop
  --------------------------------------------------------------------*/
  for (it = 0; it < maxit; it++) {
    
    if (precon) 
      precon(r, z, params);
    //    else 
    //  BLAS_CCOPY(&n, r, &ONE, z, &ONE);
   
    rhoprev = rho;
    tempc = wrap_zsum_cdot(&n, r, &ONE, z, &ONE, params); /* Sum in double */
    //tempc = wrap_cdot(&n, r, &ONE, z, &ONE, params);      /* Sum in single */
    rho = tempc.r;
    reshist[it] = sqrt(rho);
    if (plvl >= 2) printf(" Linsys res( %d ): %g\n",*iter+it,reshist[it]);

    /* Convergence test */
    if (reshist[it] < tolb) {  
       *flag = 0;
       break;  /* break do not return */
    }
    /* Restart test */
    if (reshist[it] < restartTol*reshist[0] ) {  
       *flag = 3;
       break;  /* break do not return */
    }

    if (it == 0)
      BLAS_CCOPY(&n, z, &ONE, p, &ONE);
    else {
      betaprev = beta;
      beta = rho / rhoprev;
      if (beta == 0.0) {
	*flag = 2;
	break;
      }
      for (i = 0; i < n; i++)	{
	p[i].r = z[i].r + beta * p[i].r;
	p[i].i = z[i].i + beta * p[i].i;
      }
    }

    /*----- eigCG specific code -------------------------------------------*/
    /* Remember Ap from previous iteration to be used at restart */
    
    if (nev > 0 && v_size == v_max) 
      BLAS_CCOPY(&n, Ap, &ONE, Ap_prev, &ONE);
    /*---------------------------------------------------------------------*/

    matvec(p, Ap, params);

    /*----- eigCG specific code -------------------------------------------*/
    			/* Note: I must also use the easy way to 
			 * obtain ev residual norms */
    if (nev > 0) {
       /* record the diagonal vAv for the previous std::vector */
       if (it > 0) {
	  H[(v_size-1)*v_max+v_size-1].r = 1.0/alpha + betaprev/alphaprev;
	  H[(v_size-1)*v_max+v_size-1].i = 0.0;
       }
       
       /* Restarting V */
       if (v_size == v_max) {
	  /* Solve (v_max) and (v_max-1) eigenproblems */
	  int info, allelems = v_max*v_max;
	  tempi = v_max;
	  BLAS_ZCOPY(&allelems, H, &ONE, Hevecs, &ONE);
	  BLAS_ZHEEV(&cV,&cU,&tempi,Hevecs,&v_max,Hevals,
			  zwork,&lwork,rwork,&info);
	  if (plvl < 4) {
	     /* Compute the Lanczos residual norm estimates */
	     tempd = sqrt(beta)/alpha;
	     for (i=1;i<=nev;i++) 
	         rnorms[i-1] = tempd*z_abs_primme(Hevecs[i*v_max-1]);
	  }

	  tempi = v_max-1;
	  BLAS_ZCOPY(&allelems, H, &ONE, Hevecsold, &ONE);
	  BLAS_ZHEEV(&cV,&cU,&tempi,Hevecsold,&v_max,Hevalsold,
			  zwork,&lwork,rwork,&info);
	  /* fill 0s in vmax-th elem of oldevecs to match Hevecs */
          for(i=1; i <= v_max ; i++)
	     {Hevecsold[i*v_max-1].r = 0.0 ; Hevecsold[i*v_max-1].i = 0.0;}
	  /* Attach the first nev oldevecs at the end of the nev latest ones */
	  tempi = nev*v_max;
	  BLAS_ZCOPY(&tempi,Hevecsold,&ONE,&Hevecs[tempi],&ONE);
          /* Orthogonalize the 2*nev (new+old) vectors Hevecs=QR */
	  v_size = 2*nev; 
	  BLAS_ZGEQRF(&v_max,&v_size,Hevecs,&v_max,TAU,zwork,&lwork,&info) ; 
	  /* use as a temp space Hevecsold = Q^THQ */
	  BLAS_ZCOPY(&allelems,H,&ONE,Hevecsold,&ONE); 
	  BLAS_ZUNMQR(&cR,&cN,&v_max,&v_max,&v_size,Hevecs,&v_max,
		     TAU,Hevecsold,&v_max,zwork,&lwork,&info);
	  BLAS_ZUNMQR(&cL,&cC,&v_max,&v_size,&v_size,Hevecs,&v_max,
		     TAU,Hevecsold,&v_max,zwork,&lwork,&info);
	  /* solve the small Hevecsold v_size x v_size eigenproblem */
	  BLAS_ZHEEV(&cV,&cU,&v_size,Hevecsold,&v_max,Hevals,
			  zwork,&lwork,rwork,&info);

	  /* zero out unused part of eigenectors in Hevecsold */
	  tempi = 0;
	  for(i = 0; i < v_size; i++ ) {
	     for(j = v_size; j < v_max; j++)
	        {Hevecsold[tempi + j].r=0.0; Hevecsold[tempi + j].i=0.0;}
	     tempi += v_max;
	  }
	  /* Compute the Hevecsold = Hevecs*Hevecsold */
	  BLAS_ZUNMQR(&cL,&cN,&v_max,&v_size,&v_size,Hevecs,&v_max,
		     TAU,Hevecsold,&v_max,zwork,&lwork,&info);
	  /* Copy Hevecsold into single precision */
	  for (i=0; i < v_size*v_max; i++)
	      {Coef[i].r = Hevecsold[i].r; Coef[i].i = Hevecsold[i].i;}
	  /* Restart V = V(n,v_max)*Coef(v_max,v_size) */
	  Crestart_X(V, n, Coef, n, v_max, v_size, tmp_work, tmpsize);
	  /* Restart H = diag(Hevals) plus a column and a row */
	  for (i = 0; i < allelems; i++ )  {H[i].r = 0.0; H[i].i=0.0;}
    	  for (i = 0; i < v_size; i++) H[i*(v_max+1)].r = Hevals[i];

	  if (plvl > 3) { 
	     /* Compute actual residual norms -- Expensive!*/
	     for (i=0;i<nev;i++) {
                 computeResNorm(&V[i*n], &evals[i], &rnorms[i], tmp_work, 
				 n, matvec, params);
	     }
	  }
	  if (plvl >= 3) {       /* Detailed output */
	     /* Reporting residual norm (estimates) */
	     for (i=0;i<nev;i++) 
		 printf(" Eigenvalue(%d)=%le, eigenresidual norm = %le\n",
		        i,Hevals[i],rnorms[i]);
	  }
	  
          /* The next (preconditioned) residual to be added (v = z/sqrt(rho)) 
     	   * needs the (nev+1)-th column and row, through V(:,1:vs)'*A*v. 
	   * Instead of a matvec, we use the Ap and Ap_prev to obtain this:
	   * V(:,1:vs)'*A*V(:,vs+1) = V(:,1:vs)'*A*z/sqrt(rho) = 
	   * V'(A(p-beta*p_prev))/sqrt(rho) = V'(Ap - beta*Ap_prev)/sqrt(rho)*/
    	  for (i = 0; i < n; i++) 
	      {Ap_prev[i].r = Ap[i].r-beta*Ap_prev[i].r;
	       Ap_prev[i].i = Ap[i].i-beta*Ap_prev[i].i; }
	  /* The wrap_cgemv is BLAS + a globalsum. Argum.tmp_work not in BLAS */
	  wrap_cgemv(&cC,&n,&v_size,&tpone,V,&n,Ap_prev,&ONE,
		     &tzero,Coef,&ONE,tmp_work,params);
	  /* Copy Coef into double precision H. reshist[it] holds sqrt(rho) */
	  /* Note this is not needed if instead of cgemv we use double 
	   * precision inner products (e.g., from QDP++) */
	  tempi = v_size*v_max;
	  for (i=0; i<v_size; i++) {
	      tempz.r = Coef[i].r/reshist[it];
	      tempz.i = Coef[i].i/reshist[it];
	      H[tempi + i].r = tempz.r;
	      H[tempi + i].i = tempz.i;
	      d_cnjg_primme(&H[i*v_max + v_size], &H[tempi+i]);
	  }
       } /* end of if v_size == v_max */
       else {
	  /* update (vs+1,vs),(vs,vs+1) elements of tridigonal which are real*/
          if ( it > 0) {
             H[(v_size-1)*v_max + v_size].r = -sqrt(beta)/alpha;
             H[(v_size-1)*v_max + v_size].i = 0.0;
	     H[v_size*v_max + v_size-1].r = H[(v_size-1)*v_max + v_size].r;
	     H[v_size*v_max + v_size-1].i = H[(v_size-1)*v_max + v_size].i;
	  }
	  
       } /* of else */

       /* Augment V with the current CG residual r normalized by sqrt(rho) */
       tempc.r = 1.0/reshist[it];
       tempi = v_size*n;
       for (i=0; i<n; i++) {
	   V[tempi + i].r = z[i].r*tempc.r;
	   V[tempi + i].i = z[i].i*tempc.r;
       }
       v_size++;

    } /* end of if nev >0 , ie., the eigCG specific code */
    /*---------------------------------------------------------------------*/

    /* pAp = p' * Ap */
    tempc = wrap_zsum_cdot(&n, p, &ONE, Ap, &ONE, params); /*Sum in double */
    //tempc = wrap_cdot(&n, p, &ONE, Ap, &ONE, params);  /* Sum in single */
    pAp = tempc.r;
    if (pAp == 0.0) {
      *flag = 2;
      break;
    } 

    alphaprev = alpha;
    alpha = rho / pAp;
#define USE_DOUBLE_PREC_SUM_NOT_CAXPY
#ifdef USE_DOUBLE_PREC_SUM_NOT_CAXPY
    /* The following addition is performed in double precision */
    for (i=0;i<n;i++) {
      x[i].r = x[i].r+alpha*p[i].r;   /* x = x + alpha * Ap */
      x[i].i = x[i].i+alpha*p[i].i;
      r[i].r = r[i].r-alpha*Ap[i].r;  /* r = r - alpha * Ap */
      r[i].i = r[i].i-alpha*Ap[i].i;
    }
#else
    /* Use BLAS in signle precision */
    tempc.r = alpha; tempc.i = 0.0;
    BLAS_CAXPY(&n, &tempc, p, &ONE, x, &ONE);  /* x = x + alpha * Ap */
    tempc.r = -alpha;
    BLAS_CAXPY(&n, &tempc, Ap, &ONE, r, &ONE); /* r = r - alpha * Ap */
#endif
    
//printf("%d beta, alpha, rho, pAp %le %le %le %le\n",it,beta,alpha,rho,pAp);
  } /* for it = 0 : maxit-1 */
  
  *iter = *iter + it+1; /* record the number of CG iterations plus any older */

  if (plvl) displayInfo(tol,maxit,*flag,*iter-1,reshist[it]);

  if (nev > 0) {
     /* Restart V, compute and return most recent nev Ritz pairs */
     v_size--;   /* last std::vector not projected yet */
     computeFinalEvecs(V, n, Hevals, H, v_max, v_size, nev, Coef, 
		     tmp_work, tmpsize, zwork, lwork, rwork);
     /* Compute actual residual norms -- Uses matvecs */
     for (i=0;i<nev;i++) 
         computeResNorm(&V[i*n], &evals[i], &rnorms[i], 
	             tmp_work, n, matvec, params);
  }

  free(H);
  free(Hevecs);
  free(Hevecsold);
  free(Hevals);
  free(Hevalsold);
  free(TAU);
  free(zwork);
  free(rwork);
  free(Coef);
} 
/* end of EIGPCG ************************************************************/

/* Some other functions */

/* ResNorm *******************************************************************/
/* Given a std::vector x, it normalizes it and computes and returns the
 * its residual std::vector in Res, its norm in rnorm, and its Ritz value in lambda.
 * Requires matvec() and whatever params are needed in matvec */
void computeResNorm
    (Complex_C *x, float *lambda, float *rnorm, Complex_C *Res, int n, 
    void (*matvec)(void *, void *, void *), void *params)
{
   int i, ONE = 1;
   Complex_C tempc;
   float normx2;

   /* Norm of x squared */
   tempc = wrap_zsum_cdot(&n, x, &ONE, x, &ONE, params); /* Sum in double */
   normx2 = tempc.r; 

   /* compute Ax */
   matvec(x,Res,params);

   /* lambda = x'Ax/x'x */
   tempc = wrap_zsum_cdot(&n, x, &ONE, Res, &ONE, params); /* Sum in double */
   *lambda = tempc.r/normx2;

   /* res = w - lambda v */
   for (i=0; i<n; i++) {
        Res[i].r = Res[i].r-(*lambda)*x[i].r;
        Res[i].i = Res[i].i-(*lambda)*x[i].i;
   }
   tempc = wrap_cdot(&n, Res, &ONE, Res, &ONE, params); /* Sum in single */
   *rnorm = sqrt(tempc.r/normx2);
}
/* end ResNorm ***************************************************************/

/* computeFinalEvecs *********************************************************/
/* Given V a basis of v_size and a v_size x v_size matrix H = V'*A*V, 
 * compute the nev Ritz vectors and Ritz values */
void computeFinalEvecs(Complex_C *V, int n, double *Hevals, Complex_Z *H, 
     int v_max, int v_size, int nev, Complex_C *Coef, Complex_C *tmp_work, 
     int tmpsize, Complex_Z *zwork, int lwork, double *rwork)
{
  char cV = 'V'; char cU = 'U'; 
  int i,info;
   
  /* Solve the vsize x vsize eigenproblem H */
  BLAS_ZHEEV(&cV,&cU,&v_size,H,&v_max,Hevals,zwork,&lwork,rwork,&info);

  /* Copy H small eigenvectors into single precision */
  for (i=0; i < nev*v_max; i++)
      {Coef[i].r = H[i].r; Coef[i].i = H[i].i;}
  /* Restart V = V(n,v_size)*Coef(v_max,v_size) */
  Crestart_X(V, n, Coef, n, v_max, nev, tmp_work, tmpsize);

}
/* end computeFinalEvecs ****************************************************/

/* ITERMSG - print information about iteration */
static void displayInfo(float tol,
		    int maxit,
		    int flag,
		    int iter,
		    float resnorm) {
  if (flag != 0) {
    printf("PCG stopped at iteration %d with flag %d. ", iter, flag);
  }
  
  switch(flag) {
  case 0:
    if (iter == 0)
      printf("The initial guess has relative residual %0.2g which is within\nthe desired tolerance %0.2g\n", resnorm, tol);
    else
      printf("PCG converged at iteration %d to a solution with residual norm %0.2g", iter, resnorm);
    break;
  case 1:
    printf("\nbecause the maximum number of iterations was reached.");
    break;
  case 2:
    printf("\nbecause a scalar quantity became too small.");
    break;
  }
  
  if (flag != 0)
    printf("\nThe iterate returned at iteration %d has residual norm %0.2g",iter,resnorm);

  printf("\n");
}
