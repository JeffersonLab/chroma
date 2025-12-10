// $Id: conjugate.h,v 1.1 2009-10-21 20:50:56 kostas Exp $
/*------------------------------------------------------------------------
  Determine if two numbers are complex conjugate pairs  
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/
#ifndef CONJUGATE_H_
#define CONJUGATE_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "qdp-lapack_Complex.h"

int Zconjugate(Complex_Z a, Complex_Z b, double epsi);
/* Determine if two complex numbers a and b are conjugate in an eps way. If the two numbers
   are conjugate, the return value is 1, otherwise, the return value is zero.
-------------------------------------------------------------------------------------------*/

int Cconjugate(Complex_C a, Complex_C b, float epsi);
/* same as Zconjugate for single precision */

#endif
