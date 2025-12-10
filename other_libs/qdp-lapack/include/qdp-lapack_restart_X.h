/*******************************************************************************
 * Subroutine restart_X - This subroutine computes X*hVecs and places 
 *    the result in X.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ldx 		Leading dimension of X
 *
 * nLocal       Number of rows of V assigned to the node
 *
 * basisSize    Current size of the basis V
 *
 * restartSize  Number of Ritz vectors V/W will be restarted with 
 *
 * rwork        Work array that must be at least of size restartSize
 *
 * rworkSize    The size availble in rwork. Matrix multiply blocks of X with 
 * 		hVecs, producing blocks of the new X of size
 * 		(AvailRows * restartSize) = rworkSize
 * 		Therefore rworkSize must be at least restartSize.
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * X      Holds either V or W before and after restarting
 *
 * hVecs  The eigenvectors of V'*A*V before and after restarting
 *
 ******************************************************************************/

#ifndef RESTART_X_H_
#define RESTART_X_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical_private.h"

void Crestart_X(Complex_C *X, int ldx, Complex_C *hVecs, int nLocal, 
                int basisSize, int restartSize, Complex_C *rwork, int rworkSize);



void Zrestart_X(Complex_Z *X, int ldx, Complex_Z *hVecs, int nLocal, 
                int basisSize, int restartSize, Complex_Z *rwork, int rworkSize);

#endif
