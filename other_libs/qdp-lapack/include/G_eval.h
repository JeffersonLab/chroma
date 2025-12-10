// $Id: G_eval.h,v 1.1 2009-10-21 20:50:56 kostas Exp $
/*------------------------------------------------------------------------
  Computing and sorting eignvalues and eignvectors for a general complex
  matrix.
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/

#ifndef G_EVAL_H_
#define G_EVAL_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical_private.h"
#include "sort.h"


void ZG_eval(Complex_Z *A, int N, int LDA, Complex_Z *W, char SRT_OPT, double epsi, 
             Complex_Z *VL, int LDVL, Complex_Z *VR, int LDVR, int *info);
/* Computes and sorts eigenvalues and eigenvectors of a complex double general matrix 
   which is non-defective using the LAPACK subroutine ZGEEV.

   A: is the input matrix of dimension A(LDA,N) and is stored column by cloumn 
      in a single dimension array of length LDA*N. Only the first N elements of
      each column is active. On exit, A will be destroyed.

   N: Dimension of the input matrix.
 LDA: Leading dimension of A. LDA >= N.

   W: N dimensional array of the computed eigenvalues sorted according to SRT_OPT
      as will be explained.
SRT_OPT:  Determines how to sort the eigenvalues. If SRT_OPT='R' then the eigenvalues
          are sorted in ascending order according to the smallest real part. If
          SRT_OPT='M', the eigenvalues will be sorted in ascending order according
          to having the smallest absolute value (magnitude) and such that for complex
          conjugate pairs of eigenvalues, the eigenvalue with positive imaginary part
          is listed first.The conjugacy is determined based on the parameter 'epsi'.

epsi  :   threshold used to determine if two complex numbers x and y are conjugate pairs.

          if( imag(x)*imag(y) < 0  and
              abs(imag(x)+imag(y))/abs(x+y) < epsi and
              abs(real(x)-real(y)) / abs(x+y) < epsi )

           then x and y are considered to be conjugate, otherwise they are not.


  VL: Matrix of dimesnion (LDVL,N) storing the N left eigenvectors in a one-dimensional
      array of length LDVL*N column by column where each has length LDVL such that only
      the first N elements of each column is active. The left eigenvectors satisfy:
               A^H * VL(:,k) = conjugate(W(k)) * VL(:,k)

LDVL: Leading dimension for the left eigenvectors VL such that LDVL >= N.
    
VR  : Matrix of dimension (LDVR,N) storing the N right eigenvectors in a one-dimensional
      array of length LDVR*N column by column where each column has length LDVR such that
      only the first N elements of each column is active. The right eigenvectors satisfy:
              A * VR(:,k) = W(k) * VR(:,k)

LDVR: Leading dimension of the right eigenvectors VR such that LDVR >= N.

info: Gives information about the exit status of the function based on the exit
      status of the LAPACK subroutine ZGEEV:

      info=0 successful exit
      info<0: if info=-i, the i-th argument has illegal value.
      info>0: if info=i, the QR algorithm failed to compute all eigenvalues, and no
              eigenvectors have been computed; elements and i+1:N of W contain eigenvalues
              which have converged.
=====================================

NOTE: left and right eigenvectors are bi-orthogonal with VL(:,k)^H * VR(:,j) = delta(k,j)
      where delta(k,j)=1 if k=j and is zero otherwise.
      VR and VL are normalized in a symmetrical way.

------------------------------------------------------------------------------------------*/

void CG_eval(Complex_C *A, int N, int LDA, Complex_C *W, char SRT_OPT, float epsi,  
             Complex_C *VL, int LDVL, Complex_C *VR, int LDVR, int *info);

/* Same as ZG_eval for single precision */


void ZG_eval_original(Complex_Z *A, int N, int LDA, Complex_Z *W, char SRT_OPT, double epsi, 
             Complex_Z *VL, int LDVL, Complex_Z *VR, int LDVR, int *info);




#endif
