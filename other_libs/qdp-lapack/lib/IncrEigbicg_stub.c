/*******************************************************************************
 * Function IncrEigpcg -- Incremental eigpcg. 
 *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "qdp-lapack_IncrEigpcg.h"


/******************************************************************************/

void IncrEigbicg_C( int n, int lde,int nrhs, Complex_C *X, Complex_C *B, int *ncurEvals,         	
                    int ldh, Complex_C *evecsl, Complex_C *evecsr, Complex_C *evals, 		
	            Complex_C *H,void (*matvec) (void *, void *, void *),  
	            void (*mathvec)(void *, void *, void *),void *params, float *AnormEst, 
	            Complex_C *work, Complex_C *VL,int ldvl,Complex_C *VR, int ldvr,        
	            Complex_C *ework, int esize,float tol,float *restartTol, 	
	            int maxit, char SRT_OPT, float epsi, int ConvTestOpt, int plvl,int nev,
                    int v_max,FILE *outputFile)      
{
  fprintf(stderr, "%s: not implemented in this stub", __func__);
  exit(1);
}


void IncrEigbicg_Z(  int n, int lde,int nrhs, Complex_Z *X, Complex_Z *B, int *ncurEvals,         	
                     int ldh, Complex_Z *evecsl, Complex_Z *evecsr, Complex_Z *evals, 		
	             Complex_Z *H, void (*matvec) (void *, void *, void *),  
	             void (*mathvec)(void *, void *, void *), void *params, double *AnormEst, 
	             Complex_Z *work, Complex_Z *VL, int ldvl, Complex_Z *VR, int ldvr,        
	             Complex_Z *ework, int esize, double tol, double *restartTol, 	
	             int maxit, char SRT_OPT, double epsi, int ConvTestOpt, int plvl, int nev,
                     int v_max,FILE *outputFile)
{
  fprintf(stderr, "%s: not implemented in this stub", __func__);
  exit(1);
}
/******************************************************************************/
/* END OF IncrEigpcg */
/******************************************************************************/
