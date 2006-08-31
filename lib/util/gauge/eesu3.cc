// -*- C++ -*-
// $Id: eesu3.cc,v 3.1 2006-08-31 01:31:10 bjoo Exp $
/*! \file
 *  \brief Exactly exponentiate a SU(3) lie algebra element
 */

#include "chromabase.h"
#include "util/gauge/eesu3.h"

namespace Chroma 
{

  //! Exact exponentiation of SU(3) matrix.
  /*!
   *  Input: 3x3 anti-Hermitian, traceless matrix iQ
   *  Output: exp(iQ)
   *
   *  Formula for exp(iQ) = f0 + f1*Q + f2*Q*Q is found
   *  in section III of hep-lat/0311018.
   */
  void eesu3(LatticeColorMatrix & iQ)
  {
    START_CODE( );

    LatticeColorMatrix Q = timesMinusI(iQ);

    LatticeComplex f0, f1, f2;
  
    LatticeColorMatrix QQ = Q*Q;

    LatticeReal c0    = real((1.0/3.0) * trace(Q*QQ));
    LatticeReal c1    = real((1.0/2.0) * trace(QQ));
    LatticeReal c0abs = fabs(c0);
    LatticeReal c0max = 2 * pow((c1 / 3.0), 1.5);
    LatticeReal theta = acos(c0abs/c0max);
    LatticeReal u     = sqrt(c1 / 3.0) * cos(theta / 3.0);
    LatticeReal w     = sqrt(c1) * sin(theta / 3.0);
    LatticeReal uu    = u*u;
    LatticeReal ww    = w*w;
    LatticeReal cosu  = cos(u);
    LatticeReal cosw  = cos(w);
    LatticeReal sinu  = sin(u);
    LatticeReal sinw  = sin(w);

    // exp(2iu) and exp(-iu)
    LatticeComplex exp2iu = cmplx((2*cosu*cosu - 1), 2*cosu*sinu);
    LatticeComplex expmiu = cmplx(cosu, -sinu);

    LatticeBoolean latboo_c0 = (c0      <      0);
    LatticeBoolean latboo_c1 = (c1      > 1.0e-4);
    LatticeBoolean latboo_w  = (fabs(w) >   0.05);

    LatticeReal denom = where(latboo_c1,
			      9 * uu - ww,
			      Real(1.0));

    // xi0 = xi0(w).  Expand xi0 if w is small.
    LatticeReal xi0 = where(latboo_w, 
			    sinw/w, 
			    1 - (1.0/6.0)*ww*(1-(1.0/20.0)*ww*(1-(1.0/42.0)*ww)));

    // f_i = f_i(c0, c1). Expand f_i by c1, if c1 is small.
    f0 = where(latboo_c1,
	       ((uu - ww) * exp2iu + expmiu 
		* cmplx(8*uu*cosw, 2*u*(3*uu+ww)*xi0))/denom,
	       cmplx(1-c0*c0/720, -c0/6*(1-c1/20*(1-c1/42))));

    f1 = where(latboo_c1,
	       (2*u*exp2iu - expmiu * cmplx(2*u*cosw, (ww-3*uu)*xi0))/denom,
	       cmplx(c0/24*(1.0-c1/15*(1-3*c1/112)),
		     1-c1/6*(1-c1/20*(1-c1/42))-c0*c0/5040));
  
    f2 = where(latboo_c1,
	       (exp2iu - expmiu * cmplx(cosw, 3*u*xi0))/denom,
	       0.5*cmplx(-1+c1/12*(1-c1/30*(1-c1/56))
			 +c0*c0/20160, c0/60*(1-c1/21*(1-c1/48))));
  	     
    // f_j(-c0, c1) = (-1)^j f*_j(c0, c1)
    f0 = where(latboo_c0 && latboo_c1,
	       adj(f0),
	       f0);

    f1 = where(latboo_c0 && latboo_c1,
	       -1.0*adj(f1),
	       f1);

    f2 = where(latboo_c0 && latboo_c1,
	       adj(f2),
	       f2);
  
    // evaluate f0 + f1 Q + f2 QQ (= exp(iQ)) back into iQ
    iQ = f0 + f1 * Q + f2 * QQ;

    // Relpace Q by exp(iQ)
    // Q = expiQ;


    END_CODE( );
  }

}  // end namespace Chroma


  
