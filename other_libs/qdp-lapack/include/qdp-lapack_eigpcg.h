#ifndef __EIGPCG_C_H
#define __EIGPCG_C_H

#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical.h"
#include "qdp-lapack_numerical_private.h"
#include "qdp-lapack_restart_X.h"
#include "qdp-lapack_ortho.h"

#ifdef __cplusplus
extern "C" {
#endif 

/* function prototypes */
static void displayInfo(float tol, int maxit, int flag, int iter, 
		float resnorm);

void computeResNorm
    (Complex_C *x, float *lambda, float *rnorm, Complex_C *Res, int n, 
    void (*matvec)(void *, void *, void *), void *params);

void computeFinalEvecs(Complex_C *V, int n, double *Hevals, Complex_Z *H, 
     int v_max, int v_size, int nev, Complex_C *Coef, Complex_C *tmp_work, 
     int tmpsize, Complex_Z *zwork, int lwork, double *rwork);

void eigpcg(int n, int lde, Complex_C *x, Complex_C *b, 
	 float *normb, float tol, float restartTol, int maxit,
	 int *iter, float *reshist, int *flag, int plvl, Complex_C *work,
	 void (*matvec)(void *, void *, void *),
	 void (*precon)(void *, void *, void *),
	 void *params,
	 int nev, float *evals, float *rnorms,
  	 int v_max, Complex_C *V, int esize, Complex_C *ework);

#ifdef __cplusplus
};
#endif 
#endif
