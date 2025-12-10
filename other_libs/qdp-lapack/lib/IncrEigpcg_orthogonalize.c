/* $Id: IncrEigpcg_orthogonalize.c,v 1.2 2009-10-21 20:50:57 kostas Exp $ */
/*******************************************************************************
 * Function IncrEigpcg -- Incremental eigpcg. 
 *
 ******************************************************************************/
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "qdp-lapack_IncrEigpcg.h"
#include "qdp-lapack_wtime.h"

/* DEBUG **********************************************************************/
void RayleighRitz(Complex_C *evecs, int lde, int n, int numEvecs,
     Complex_C *H, int ldh,
     FILE *outputFile, void (*matvec)(void *, void *, void *),  void *params)
{
     int Eelems = lde*numEvecs;
     int Helems = numEvecs*numEvecs;
     int lwork = 2*ldh;
     int ONE = 1;
     int i,info;
     Complex_C *V;
     Complex_C *Ht;
     Complex_C *cwork;
     float *Hevals;
     float *awork;
     if (lwork < n) lwork = n;

     float rnorm;
     char cV = 'V'; char cU = 'U'; 

     V = (Complex_C *) calloc(Eelems, sizeof(Complex_C));
     Ht = (Complex_C *) calloc(Helems, sizeof(Complex_C));
     cwork = (Complex_C *) calloc(lwork, sizeof(Complex_C));
     Hevals = (float *) calloc(numEvecs, sizeof(float));
     awork = (float *) calloc(3*ldh, sizeof(float));
	     
     BLAS_CCOPY(&Eelems, evecs, &ONE, V, &ONE);
     for (i=0;i<numEvecs;i++)
        BLAS_CCOPY(&numEvecs, &H[i*ldh], &ONE, &Ht[i*numEvecs], &ONE);

     BLAS_CHEEV(&cV,&cU,&numEvecs,Ht,&numEvecs,Hevals,cwork,&lwork,awork,&info);
     Crestart_X(V, lde, Ht, n, numEvecs, numEvecs, cwork, lwork);

     for (i=0;i<numEvecs;i++)  {
         computeResNorm(&V[i*lde], &Hevals[i], &rnorm, cwork, n, matvec,params);
         fprintf(outputFile, "RR Eval[%d]: %22.15E rnorm: %22.15E\n", 
			                i+1, Hevals[i], rnorm);
     }
	
     free(V);
     free(Ht);
     free(cwork);
     free(Hevals);
     free(awork);
}

/* end DEBUG ******************************************************************/

/******************************************************************************/
void IncrEigpcg(int n, int lde, /* n dim of matrix A, lde leading dim of vecs */
	int nrhs, 		/* The number of right hand sides to solve */
	Complex_C *X,   	/* The nrhs solutions */
	Complex_C *B,   	/* The nrhs right hand sides */
	int *ncurEvals,		/* Num of eigenpairs already in evecs,evals */
	int ldh,  		/* Max eigenpairs that can be stored in evecs */
	Complex_C *evecs, 	/* The eigenstd::vector basis augmented by eigpcg */
	float *evals, 		/* The eigenvalues augmented by eigpcg */
	Complex_C *H,		/* The ncurEvals^2 matrix: H=evecs'*A*evecs; */
	Complex_C *HU,  	/* The Cholesky factorization of the H = U'U */
				/* ---------------------------------------- */
	void (*matvec)(void *, void *, void *),  /* Matvec operator */
	void (*precon)(void *, void *, void *),  /* Precondition operator */
	void *params,  		/* parameters to be passed to operators */
				/* ---------------------------------------- */
			/* Work arrays will be allocated if null on input */
	Complex_C *work,        /* work array for CG etc (max(4lde,nev*ldh) */
	Complex_C *V, 		/* work array for eigenstd::vector search basis */
	Complex_C *ework, 	/* work array. If esize < N+2*nev allocate */
	int esize, 		/* N+2*nev <= esize <= (2*nev+1)*N */
				/* ---------------------------------------- */
	float tol, 		/* CG linear system tolerance |res| < tol*b */
	float *restartTol, 	/* Restart CG when |res|<restartTol*|b-Ax(0)| */
	float normAestimate,    /* An estimate of ||A||_2 */
	int updateRestartTol,   /* Whether to update restartTol from eresids
				 * Expensive. Requires computation of residuals
				 * =0 Never update restartTol
				 * =1 Compute all eigenresiduals and update 
				 *    when ncurEvals=ldh for the first time
				 * =2 Update based on up to 10 eres picked from 
				 *    ncurEvals, on every rhs that adds evecs 
				 * =3 Compute all eres and update on every rhs
				 *    that adds evecs--unnecessarily expensive
				 * If updateRestartTol>0 Cholesky is not used */
	int maxit, 		/* Maximum number of CG iterations per system */
	int plvl, 		/* Printing level [0..5] */
	int nev, 		/* Number of eigenvalues to target in each CG */
	int v_max,		/* Maximum number of n-dim vectors in V */
	FILE *outputFile)       /* File to write reports. Could be STDERR */
{

  /* Timing vars */
  double wt1,wt2,wE,wI;

  /* Pointers */
  Complex_C *oldEwork = ework;
  Complex_C *x, *b;
  float *rnorms, *reshist;
 
  /* Variables */
  float machEps = 2e-8;
  float lambda, normb;
  Complex_C tempc;
  int i,j, ONE = 1;
  int zs, cs, ds, tmpsize;
  int freeV = 0;
  int freeWork = 0;
  int freeEwork = 0;
  int numIts, flag, nAdded, nev_used, iters_used;
  int maxit_remain, normsComputed, step;
  int Cholesky = 1;

  Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
  char cR = 'R'; char cL = 'L'; char cN ='N'; 
  char cV = 'V'; char cU = 'U'; char cC ='C';

  if (updateRestartTol) Cholesky = 0;

  /* ------------------------------------------------------------------- */
/* printf("--------------------------------------------\n");
printf("n, lde, nrhs,*ncurEvals,ldh,esize,*restartTol,normAestimate, updateRestartTol, nev, v_max, %d %d %d %d %d %d %g %g %d %d %d \n",
	n, lde, nrhs,*ncurEvals,ldh,esize,*restartTol,
        normAestimate, updateRestartTol, nev, v_max);
*/

  /* ------------------------------------------------------------------- */
  /* Work allocations */
  /* ------------------------------------------------------------------- */
  zs = sizeof(Complex_Z); 
  cs = sizeof(Complex_C); 
  ds = sizeof(double);
  if (work == NULL) {
     tmpsize = 4*lde;
     if (tmpsize < nev*ldh)
	tmpsize = nev*ldh;
     if ((work = calloc(tmpsize, cs)) == NULL) 
        fprintf(stderr, "ERROR IncrEigpcg could not allocate work\n");
     freeWork = 1;
  }
  if (V == NULL && nev > 0) {
     if ((V = calloc(v_max*n, cs)) == NULL) 
        fprintf(stderr, "ERROR IncrEigpcg could not allocate V\n");
     freeV = 1;
  }
  if (ework == NULL && nev > 0) {
     if (esize < n+2*nev)
	esize = 3*n;      /* Defaulting to 3*nev */
     if ((ework = calloc(esize, cs)) == NULL) 
        fprintf(stderr, "ERROR IncrEigpcg could not allocate ework\n");
     freeEwork = 1;
  }
  else { /* Not null but check if enough space */
     if (esize < n+2*nev && nev > 0) { /* Realloc ework but remember old pntr */
        freeEwork = 1;
	esize = 3*n;      /* Defaulting to 3*nev and allocating */
        if ((ework = calloc(esize, cs)) == NULL) 
           fprintf(stderr, "ERROR IncrEigpcg could not allocate ework\n");
     }
  }
  if (  (rnorms = calloc(ldh, sizeof(float)) ) == NULL )
     fprintf(stderr, "ERROR IncrEigpcg could not allocate rnorms\n");
  if (  (reshist = calloc(maxit, sizeof(float)) ) == NULL)
     fprintf(stderr, "ERROR IncrEigpcg could not allocate reshist\n");

  /* ------------------------------------------------------------------- */
  /* end Work allocations */
  /* ------------------------------------------------------------------- */

  if (*ncurEvals > 0) {
    
    /* V = A*evecs(new) */
    for (i=0; i< (*ncurEvals); i++) {  
      matvec(&evecs[i*lde], &V[0], params);

         /* Hnew = evecs(old+new)'*V(new) = evecs(old+new)'*A*evecs(new) */
         nAdded = 1;
         tmpsize = i+nAdded;
         wrap_cgemm(&cC, &cN, &tmpsize, &nAdded, &n, &tpone, evecs, &lde, 
	         V, &n, &tzero, &H[i*ldh], &ldh, work, params);
    }

    /* Copy H into HU */
    tmpsize = ldh*ldh;
    BLAS_CCOPY(&tmpsize, H, &ONE, HU, &ONE);
   
    if (Cholesky)  /* Cholesky factorize H = HU'*HU */
        BLAS_CPOTRF(&cU, ncurEvals, HU, &ldh, &flag);
    else {         /* Eigendecompose: H = HU Lambda HU'. Don't form V(HU)*/
       i = 4*lde;
       BLAS_CHEEV(&cV,&cU,ncurEvals,HU,&ldh,evals,work,&i,(float *)ework,&flag);
    }

 } /* if nev>0 */


  /* ------------------------------------------------------------------- */
  /* Solving one by one the nrhs systems with incremental init-eigpcg    */
  /* ------------------------------------------------------------------- */
  for (j=0; j<nrhs; j++) {

      /* The j-th system */
      x = &X[j*lde];
      b = &B[j*lde]; 

      if (plvl) fprintf(outputFile, "System %d\n", j);

      /* --------------------------------------------------------------- */
      /* WHILE: Call eigCG until this rhs converges. eigCG must be called 
       * repeatedly if the deflation in initCG needs to restart it often */
      /* --------------------------------------------------------------- */
      wE = 0.0; wI = 0.0;     /* Start accumulator timers */
      flag = -1;    	      /* First time through. Run eigCG regularly */
      maxit_remain = maxit;   /* Initialize Max and current # of iters   */
      numIts = 0;
      while (flag == -1 || flag == 3) {

         if (*ncurEvals > 0) {
            /* --------------------------------------------------------- */
            /* Perform init-CG with evecs vectors
             * xinit = xinit + evecs*Hinv*evec'*(b-Ax0) 		 */
            /* --------------------------------------------------------- */
            wt1 = primme_get_wtime(); 
   
            tempc = wrap_cdot(&n, x, &ONE, x, &ONE, params); /* norm ||x||^2 */
            if (tempc.r > 0.0) { /* If initial guess */
               matvec(x, work, params);  /* work(1:n) used for residual b-Ax */
               for (i = 0; i < n; i ++) {
                   work[i].r = b[i].r - work[i].r;
                   work[i].i = b[i].i - work[i].i;
               }
            }
            else 
  	       BLAS_CCOPY(&n, b, &ONE, work, &ONE);
	    
            /* work(n:n+curEvals) = evecs(n x ncurEvals)'*work(1:n) */
            wrap_cgemv(&cC,&n,ncurEvals,&tpone,evecs,&lde,work,&ONE,
	           &tzero,&work[n],&ONE,&work[2*n],params);
   
	    if (Cholesky) 
            /* HU is the upper triangular factor of the Cholesky */
            /* work(n:n+ncurEvals)=H^(-1)*work(n:n+ncurEvals) =(evecs'*work) */
	       BLAS_CPOTRS(&cU,ncurEvals,&ONE,HU,&ldh,&work[n],ncurEvals,&flag);
	    else {
            /* HU is the eigenvectors of the matrix H */
            /* work(n:n+ncurEvals) = HU*Lambda^(-1)*HU'*work(n:n+ncurEvals) */
               BLAS_CGEMV(&cC,ncurEvals,ncurEvals,&tpone,HU,&ldh,&work[n],&ONE,
		       &tzero,work,&ONE);
	       for (i = 0; i<*ncurEvals ; i++) {
	           work[i].r = work[i].r/evals[i];
	           work[i].i = work[i].i/evals[i];
	       }
               BLAS_CGEMV(&cN,ncurEvals,ncurEvals,&tpone,HU,&ldh,work,&ONE,
		       &tzero,&work[n],&ONE);
	    }
            /* work[1:n] = evecs*work(n:n+ncurEvals) */
            BLAS_CGEMV(&cN,&n,ncurEvals,&tpone,evecs,&lde,&work[n],&ONE,
		       &tzero,work,&ONE);
            /* x = x + work the new initial guess */
	    BLAS_CAXPY(&n, &tpone, work, &ONE, x, &ONE);
   
            wt2 = primme_get_wtime();
	    wI = wI + wt2-wt1;
         }
         /* end of init-CG with evecs vectors                            */
         /* ------------------------------------------------------------ */
   
         /* ------------------------------------------------------------ */
	 /* Adjust nev for eigcg according to available ldh/restart      */
         /* ------------------------------------------------------------ */
	 if (flag == 3) { /* restart with the same rhs, set nev_used = 0 */
	    nev_used = 0;
	    /* if convergence seems before next restart do not restart again */
 	    if (reshist[numIts-1]*(*restartTol) < tol*normb) 
		*restartTol = 0.0;
	 }
	 else {         /* First time through this rhs. Find nev evecs */
                        /* limited by the ldh evecs we can store in total */
            if (ldh-*ncurEvals < nev)
	       nev = ldh - *ncurEvals;
	    nev_used = nev;
	 }
         /* ------------------------------------------------------------ */
         /* Solve Ax = b with x initial guess                            */
         /* ------------------------------------------------------------ */
         wt1 = primme_get_wtime(); 

	 eigpcg(n, lde, x, b, &normb, tol, *restartTol, 
		maxit_remain, &numIts, &reshist[numIts], &flag, plvl,
		work, matvec, precon, params,
                nev_used, &evals[*ncurEvals], &rnorms[*ncurEvals], 
		v_max, V, esize, ework);

         wt2 = primme_get_wtime();
	 wE = wE + wt2-wt1;

	 /* if flag == 3 update the remain max number of iterations */
	 maxit_remain = maxit - numIts;
      } 
      /* end of while(this rhs !converged) call eigCG loop            */
      /* ------------------------------------------------------------ */

      /* ---------- */
      /* Reporting  */
      /* ---------- */

      if (plvl) {
         fprintf(outputFile, "For this rhs:\n");
         fprintf(outputFile, "Total initCG Wallclock : %-f\n", wI);
         fprintf(outputFile, "Total eigpcg Wallclock : %-f\n", wE);
         fprintf(outputFile, "Iterations: %-d\n", numIts); 
         fprintf(outputFile, "Actual Resid of LinSys  : %e\n", reshist[numIts-1]);
	 if (plvl > 1) 
            for (i=0; i < nev; i++) 
               fprintf(outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", 
	                i+1, evals[*ncurEvals+i], rnorms[*ncurEvals+i]); 
   
         if (flag != 0) {
            fprintf(outputFile, 
               "Error: eigpcg returned with nonzero exit status\n");
            return;
         }
      }
      /* ------------------------------------------------------------------- */

      /* ------------------------------------------------------------------- */
      /* Update the evecs and the factorization of evecs'*A*evecs            */
      /* ------------------------------------------------------------------- */
      wt1 = primme_get_wtime(); 

      if (nev > 0) {
         /* Append new Ritz pairs to evecs */
         for (i=0; i<nev; i++)
	     BLAS_CCOPY(&n, &V[i*n], &ONE, &evecs[(*ncurEvals+i)*lde], &ONE);
   
         /* Orthogonalize the new Ritz vectors. This is necessary only if
	  * we eigen-decompose H. With Cholesky it only reduces conditioning */
	 if (Cholesky)
            nAdded = ortho_new_vectors(evecs, lde, n, *ncurEvals, 2,
                &evecs[(*ncurEvals)*lde], nev, machEps, ework, work, params);
	 else
            nAdded = ortho_new_vectors(evecs, lde, n, *ncurEvals, 3,
                &evecs[(*ncurEvals)*lde], nev, machEps, ework, work, params);

         /* V = A*evecs(new) */
         for (i=0; i<nAdded; i++) 
	     matvec(&evecs[(*ncurEvals+i)*lde], &V[i*n], params);
   
         /* Hnew = evecs(old+new)'*V(new) = evecs(old+new)'*A*evecs(new) */
         tmpsize = (*ncurEvals)+nAdded;
         wrap_cgemm(&cC, &cN, &tmpsize, &nAdded, &n, &tpone, evecs, &lde, 
	         V, &n, &tzero, &H[(*ncurEvals)*ldh], &ldh, work, params);
   
         *ncurEvals = *ncurEvals + nAdded;
   
         /* Copy H into HU */
         for (i=0; i<*ncurEvals; i++)
     	    BLAS_CCOPY(ncurEvals, &H[i*ldh], &ONE, &HU[i*ldh], &ONE);
   
	 if (Cholesky)  /* Cholesky factorize H = HU'*HU */
            BLAS_CPOTRF(&cU, ncurEvals, HU, &ldh, &flag);
	 else {         /* Eigendecompose: H = HU Lambda HU'. Don't form V(HU)*/
	    i = 4*lde;
            BLAS_CHEEV(&cV,&cU,ncurEvals,HU,&ldh,evals,work,&i,(float *)ework,
		       &flag);
	 }

         /* ------------------------------------------------------------- */
	 /* Computing residual norms to update the new restartTol         */
	 if (updateRestartTol) {
 	    step = *ncurEvals+1;  /* so that the loop does not run */
	    switch (updateRestartTol) {
	    case 1: /* Compute all resnorms at the end of Incrementaleigcg */
		    /* Since nev>0 evecs were just added. nev=0 from now on */
	       if (*ncurEvals == ldh) step = 1; 
	       break;
	    case 2: /* Compute all resnorms on any step that adds evecs nev>0 */
	       step = 1; break;
	    case 3: /* Compute a max of 10 resnorms (with step at least nev) */
	       step = (int) (*ncurEvals/10);
	       if (step < nev) step = nev;
	       break;
	    } /* end of switch */
	    /* Note: compute norms only for evecs computed in previous rhs */
	    /* Note: do not compute the norm of the most recently added evec */
	    normsComputed = 0;
     	    for (i=step-1;i<*ncurEvals-step;i+=step)  {
		if (rnorms[i] > tol*10.0 || rnorms[i] == 0.0) {
                   BLAS_CGEMV(&cN,&n,ncurEvals,&tpone,evecs,&lde,&HU[i*ldh],
	                       &ONE,&tzero,work,&ONE);
                   computeResNorm(work, &lambda, &rnorms[i], &work[n], n, 
			       matvec, params);
		}
	        normsComputed++;
	    } /* end of for */

	    /* Update restartTol if (update every rhs) or (the last rhs)  */
	    /* Currently: arbitrary as the residual of the middle of eres */
	    /* ToDo: Optimal choice from how condition number changes     */

	    if ( (updateRestartTol != 1) || (*ncurEvals == ldh) )
	       *restartTol = rnorms[(normsComputed/2)*step]/normAestimate; 

	    if (plvl) fprintf(outputFile, "New restartTol= %f\n",*restartTol);

	 } /* end of if updateRestartTol */
         /* ------------------------------------------------------------- */

	 /* Reporting */
         wt2 = primme_get_wtime();
         if (plvl) {
            fprintf(outputFile, "Update\n");
            fprintf(outputFile, "Added %d vecs\n",nAdded);
            fprintf(outputFile, "U Wallclock : %-f\n", wt2-wt1);
	    if (plvl >= 5 ) { 
	       if (Cholesky) RayleighRitz(evecs,lde,n,*ncurEvals,H,ldh,
			       		  outputFile,matvec,params);
	       else {
     		  for (i=0;i<*ncurEvals;i++)  {
                     BLAS_CGEMV(&cN,&n,ncurEvals,&tpone,evecs,&lde,&HU[i*ldh],
			       &ONE,&tzero,work,&ONE);
         	     computeResNorm(work, &lambda, &rnorms[i],
				      &work[n], n, matvec, params);
         	     fprintf(outputFile,"RR Eval[%d]: %22.15E rnorm: %22.15E\n",
				      i+1, lambda, rnorms[i]);
		  }
	       } /* Printing if evals are known */
	    } /* plvl >=5 */
         } /* plvl > 0 */

      } /* if nev>0 */
      /* ------------------------------------------------------------------- */
      /* end of update phase */
      /* ------------------------------------------------------------------- */

   } /* end of nrhs loop */

   if (freeEwork) {
      free(ework);
      ework = oldEwork;
   }
   if (freeWork) free(work);
   if (freeV) free(V);
   free(rnorms);
   free(reshist);

   return;
}
/******************************************************************************/
/* END OF IncrEigpcg */
/******************************************************************************/
