/* $Id: restart_X.c,v 1.4 2009-10-21 20:50:57 kostas Exp $ */
/*******************************************************************************
 * Subroutine restart_X - This subroutine computes X*hVecs and places 
 *    the result in X.
 ******************************************************************************/

#include "qdp-lapack_restart_X.h"
#define min(a, b) (a < b ? a : b)


void Crestart_X(Complex_C *X, int ldx, Complex_C *hVecs, int nLocal, 
               int basisSize, int restartSize, Complex_C *rwork, int rworkSize){
   char cN = 'N';
   int ONE = 1;
   int i, k;  /* Loop variables */
   int AvailRows = min(rworkSize/restartSize, nLocal);
   Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
   i = 0;
   while (i < nLocal) {
     /* Block matrix multiply */
      BLAS_CGEMM(&cN, &cN, &AvailRows, &restartSize, &basisSize, &tpone,
         &X[i], &ldx, hVecs, &basisSize, &tzero, rwork, &AvailRows );
     /* Copy the result in the desired location of X */
      for (k=0; k < restartSize; k++) {
         BLAS_CCOPY(&AvailRows, &rwork[AvailRows*k],&ONE, &X[i+ldx*k],&ONE);
      }
      i = i+AvailRows;
      AvailRows = min(AvailRows, nLocal-i);
   }

}

/**************************************************************************************/


void Zrestart_X(Complex_Z *X, int ldx, Complex_Z *hVecs, int nLocal, 
               int basisSize, int restartSize, Complex_Z *rwork, int rworkSize){
   char cN = 'N';
   int ONE = 1;
   int i, k;  /* Loop variables */
   int AvailRows = min(rworkSize/restartSize, nLocal);
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
   i = 0;

   while (i < nLocal) {
      /* Block matrix multiply */
      BLAS_ZGEMM(&cN, &cN, &AvailRows, &restartSize, &basisSize, &tpone,
         &X[i], &ldx, hVecs, &basisSize, &tzero, rwork, &AvailRows );

      /* Copy the result in the desired location of X */
      for (k=0; k < restartSize; k++) {
         BLAS_ZCOPY(&AvailRows, &rwork[AvailRows*k],&ONE, &X[i+ldx*k],&ONE);
      }

      i = i+AvailRows;
      AvailRows = min(AvailRows, nLocal-i);
   }


}


