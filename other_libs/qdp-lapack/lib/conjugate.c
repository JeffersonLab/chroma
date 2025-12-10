/*------------------------------------------------------------------------
  Determine if two numbers are complex conjugate pairs  
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/

#include "conjugate.h"

int Zconjugate(Complex_Z a, Complex_Z b, double epsi)
{
    double ia = a.i;
    double ib = b.i;
    double ra = a.r;
    double rb = b.r;
    int flag = 0;
    Complex_Z sum;
    double sum_abs;

    sum.r = a.r+b.r;    sum.i=a.i+b.i;

    sum_abs = z_abs_primme(sum);

    if ( ia*ib < 0 ) 
      if ( (fabs(ia+ib)/sum_abs) < epsi)
        if ( (fabs(ra-rb)/sum_abs) < epsi)
	    flag = 1;
      
    return flag;
}


int Cconjugate(Complex_C a, Complex_C b, float epsi)
{
    float ia = a.i;
    float ib = b.i;
    float ra = a.r;
    float rb = b.r;
    int flag = 0;
    Complex_C sum;
    float sum_abs;

    sum.r = a.r+b.r;    sum.i=a.i+b.i;

    sum_abs = c_abs_primme(sum);

    if ( ia*ib < 0 ) 
      if ( (fabs(ia+ib)/sum_abs) < epsi)
        if ( (fabs(ra-rb)/sum_abs) < epsi)
	    flag = 1;
      
    return flag;
}



