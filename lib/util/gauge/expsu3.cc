/*#  $Id: expsu3.cc,v 3.1 2006-08-25 23:46:37 edwards Exp $ ($Date: 2006-08-25 23:46:37 $) */
/*#  EXPSU3: */
/*#  In place  a = 1 + a + (1/2)*a^2 + ...+ (1/n!)*(a)^n  n = 12 */
/*#  a must be a lattice_complex32_colour_colour  primitive */

/*#  If it does not satisfy the unitary condition, we eesu3 () to invoke */
/*#  the exact exponentiator. */
/*#  cflag decides whether we do the checking */
#include "chromabase.h" 

#include "util/gauge/expsu3.h"


namespace Chroma { 

void expsu3(LatticeColorMatrix& a, Reunitarize cflag)
{
  START_CODE();

  LatticeColorMatrix aux1;
  LatticeColorMatrix aux2;
  LatticeColorMatrix aux3;
  LatticeReal aa;
  LatticeReal bb;
  LatticeReal te1;
  LatticeReal te2;
  LatticeReal te3; 
  LatticeReal tmp;
  LatticeReal x0_real;
  LatticeReal x0_imag;
  LatticeReal x1_real;
  LatticeReal x1_imag;
  LatticeReal x2_real;
  LatticeReal x2_imag;
  LatticeComplex x0;
  LatticeComplex x1;
  LatticeComplex x2;
  LatticeBoolean bad;
  Real ftmp1;
  Real ftmp2;
  
  int numbad;
  int i;

  if ( Nc != 3 )
    QDP_error_exit("can only handle SU[3]", Nc);
  
  
  /*# aux1 = a*a */
  aux1 = a * a;
  
  /*# aa = -tr(a*a)/2 */
  aa = -Real(0.5)*real(trace(aux1));
  
  /*# aux2 = a*a*a */
  aux2 = aux1 * a;
  
  /*# bb = -Im( tr(a*a*a)/3 ) */
  bb = -imag(trace(aux2))/Real(3);
    
  ftmp1 = Real(4);
  ftmp2 = Real(396);
  te1 = aa * ftmp1;		/* multiply by 4.0 */
  tmp=ftmp2;
  te1 -= tmp;			/* subtract 396.0 */
  
  ftmp2 = Real(23760);
  te2 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te2 += tmp;				/* add 23760.0 */
  
  ftmp2 = Real(665280);
  te1 = aa * te2;			/* multiply by te1 */
  tmp=ftmp2;
  te1 -= tmp;			/* subtract 665280.0 */
  
  te2 = bb * bb;			/* multiply by bb */
  
  te3 = te1;
  te3 += te2;				/* te3 = te1 + te2 */
  
  ftmp2 = Real(479001600);
  te1 = te2 * te3;		/* multiply by te1 */
  tmp=ftmp2;
  te1 += tmp;				/* add 479001600.0 */
  
  ftmp1 = Real(2.08767569878681E-09);
  x0_real = te1 * ftmp1;		/* multiply by 2.08767569878681e-09 */
  
  
  
  ftmp1 = -Real(3);
  ftmp2 = Real(110);
  te1 = aa * ftmp1;		/* multiply by -3.0 */
  tmp=ftmp2;
  te1 += tmp;				/* add 110.0.0 */
  
  x0_imag = te1 * te2;		/* multiply by te2 */
  
  tmp=ftmp2;
  te1 = tmp;			/* te1 = 110.0 - aa */
  te1 -= aa;
  
  ftmp2 = Real(7920);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* add 7920.0 */
  
  ftmp2 = Real(332640);
  te1 = aa * te3;			/* multiply by te3 */
  tmp=ftmp2;
  te1 += tmp;				/* add 332640.0 */
  
  ftmp2 = Real(6652800);
  te3 = aa * te1;			/* multiply by te3 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 6652800.0 */
  
  te1 = te3;
  te1 += x0_imag;			/* te1 = te3 + x0_imag */
  
  te3 = bb * te1;			/* multiply by te1 */
  
  ftmp1 = Real(2.505210838544172E-08);
  x0_imag = te3 * ftmp1;		/* multiply by 2.505210838544172e-08 */
  
  
  
  ftmp1 = -Real(6);
  ftmp2 = Real(330);
  te1 = aa * ftmp1;		/* multiply by -6.0 */
  tmp=ftmp2;
  te1 += tmp;				/* add 330.0 */
  
  ftmp2 = Real(7920);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 7920.0 */
  
  x1_real = te3 * te2;		/* x1_real = te3 * te2 */
  
  ftmp2 = Real(110);
  tmp=ftmp2;
  te1 = tmp;			/* te1 = 110.0 - aa */
  te1 -= aa;
  
  ftmp2 = Real(7920);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 7920.0 */
  
  ftmp2 = Real(332640);
  te1 = aa * te3;			/* multiply by te1 */
  tmp=ftmp2;
  te1 += tmp;			 	/* add 332640.0 */
  
  ftmp2 = Real(6652800);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 6652800.0 */
  
  ftmp2 = Real(39916800);
  te1 = aa * te3;			/* multiply by te3 */
  tmp=ftmp2;
  te1 += tmp;				/* add 39916800.0 */
  
  x1_real += te1;			/* x1_real += te1 */
  
  ftmp1 = Real(2.505210838544172E-08);
  x1_real = x1_real * ftmp1;	/* multiply by 2.505210838544172e-08 */
  
  
  
  ftmp1 = -Real(4);
  ftmp2 = Real(132);
  te1 = aa * ftmp1;		/* multiply by -4.0 */
  tmp=ftmp2;
  te1 += tmp;				/* add 132.0 */
  
  x1_imag = te1 * te2;		/* x1_imag = te1 * te2 */
  
  ftmp1 = -Real(5);
  ftmp2 = Real(528);
  te1 = aa * ftmp1;		/* multiply by -5.0 */
  tmp=ftmp2;
  te1 += tmp;				/* add 528.0 */
  
  ftmp2 = Real(35640);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 35640.0 */
  
  ftmp2 = Real(1330560);
  te1 = aa * te3;			/* multiply by te1 */
  tmp=ftmp2;
  te1 += tmp;				/* add 1330560.0 */
  
  ftmp2 = Real(19958400);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 19958400.0 */
  
  te1 = x1_imag;			/* te1 = x1_imag + te3 */
  te1 += te3;
  
  ftmp1 = Real(2.08767569878681E-09);
  x1_imag = bb * te1;		/* x1_imag = bb * te1  */
  x1_imag = x1_imag * ftmp1;	/* x1_imag *= 2.08767569878681e-09 */
  
  
  
  ftmp1 = -Real(6);
  ftmp2 = Real(396);
  te1 = aa * ftmp1;		/* multiply by -6.0 */
  tmp=ftmp2;
  te1 += tmp;				/* add 396.0 */
  
  ftmp2 = Real(11880);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 11880.0 */
  
  x2_real = te3 * te2;		/* x2_real = te2 * te3 */
  
  ftmp2 = Real(132);
  te1=ftmp2;				/* te1 = 132.0 - aa */
  te1 -= aa;
  
  ftmp2 = Real(11880);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 11880.0 */
  
  ftmp2 = Real(665280);
  te1 = aa * te3;			/* multiply by te3 */
  tmp=ftmp2;
  te1 += tmp;				/* add 665280.0 */
  
  ftmp2 = Real(19958400);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 -= tmp;			/* subtract 19958400.0 */
  
  ftmp2 = Real(239500800);
  te1 = aa * te3;			/* multiply by te3 */
  tmp=ftmp2;
  te1 += tmp;				/* add 239500800.0 */
  
  ftmp1 = Real(2.08767569878681E-09);
  x2_real += te1;			/* x2_real += te1 */
  x2_real = ftmp1 * x2_real;	/* x2_real *= 2.08767569878681e-09 */
  
  
  
  ftmp1 = Real(4);
  ftmp2 = Real(330);
  te1 = ftmp1 * aa;		/* multiply by 4.0 */
  tmp=ftmp2;
  te1 -= tmp;			/* subtract 330.0 */
  
  ftmp2 = Real(15840);
  te3 = aa * te1;			/* multiply by te1 */
  tmp=ftmp2;
  te3 += tmp;				/* add 15840.0 */
  
  ftmp2 = Real(332640);
  te1 = aa * te3;			/* multiply by te3 */
  tmp=ftmp2;
  te1 -= tmp;			/* subtract 332640.0 */
  
  te3 = te1;			/* te3 = te1 + te2 */
  te3 += te2;
  
  ftmp1 = Real(2.505210838544172E-08);
  x2_imag = bb * te3;		/* x2_imag = bb * te3 */
  x2_imag = ftmp1 * x2_imag;	/* x2_imag *= 2.505210838544172e-08 */
  
  
  /*# Build lattice complex variables out of two lattice real variables. */
  x0 = cmplx(x0_real,x0_imag);
  x1 = cmplx(x1_real,x1_imag);
  x2 = cmplx(x2_real,x2_imag);
  
                
  aux2 = a * x1;
  aux2 += aux1 * x2;
  
  aux3 = 1;
  aux2 += x0 * aux3;
    
  if ( cflag == REUNITARIZE_LABEL ) {
    reunit (aux2, bad, numbad, REUNITARIZE_LABEL);
    if ( numbad > 0 ) {
	  
      QDP_error_exit("found some matrices violating unitarity", numbad);
      /*#   eesun (a, aux1, aa, bb, aux2); */
    }
  }
  a = aux2;

  END_CODE();
}


} // End namespace Chroma 
