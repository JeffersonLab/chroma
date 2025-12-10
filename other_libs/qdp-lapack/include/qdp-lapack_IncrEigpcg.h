#ifndef  __INCREIGCG_C_H
#define __INCREIGCG_C_H
#include "qdp-lapack_eigpcg.h"
#include "stdio.h"

#ifdef __cplusplus
extern "C" {
#endif 

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
	FILE *outputFile);      /* File to write reports. Could be STDERR */

#ifdef __cplusplus
};
#endif 
#endif
