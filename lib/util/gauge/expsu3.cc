// $Id: expsu3.cc,v 1.1 2004-01-22 22:52:54 edwards Exp $
/*! \file
 *  \brief Exponentiate a SU(3) matrix to 12th order
 */

#include "chromabase.h"
#include "util/gauge/expsu3.h"
#include "util/gauge/reunit.h"

using namespace QDP;

//! Exponentiate a SU(3) lie algebra element
/*!
 * \ingroup gauge
 *
 *  In place  a = 1 + a + (1/2)*a^2 + ...+ (1/n!)*(a)^n  n = 12
 *  a must be a lattice_complex32_colour_colour  primitive
 *
 *  If it does not satisfy the unitary condition, we eesu3 () to invoke
 *  the exact exponentiator.
 *  cflag decides whether we do the checking
 *
 *  Follows notes of ADK: Oct. 23, 1991, (file  exp_su3.tex) .
 *
 *  \param m        LatticeColorMatrix          (Modify)
 */

void expsu3(LatticeColorMatrix& m, int cflag)
{
  START_CODE("expsu3");
  
  if ( Nc != 3 )
    QDP_error_exit("can only handle SU[3]");
  
  // msq = m*m
  LatticeColorMatrix msq = m * m;
  
  // a = -tr(m*m)/2
  LatticeReal a = -0.5*real(trace(msq));
  
  // b = -Im( tr(m*m*m)/3 )
  Real mthird = -1.0/3.0;
  LatticeReal b = mthird*imag(trace(msq * m));
  LatticeReal bsq = b*b;
  
  // exp(M) = a2*M^2 + a1*M + a0*I
  LatticeComplex a0 = 
    cmplx(2.08767569878681e-09*(bsq*(bsq+a*(a*(4.0*a-396.0)+23760.0)-665280.0)+479001600.0),
	  2.505210838544172e-08*b*((110.0-3.0*a)*bsq+a*(a*((110.0-1.0*a)*a-7920.0)+332640.0)-6652800.0));

  LatticeComplex a1 = 
    cmplx(2.505210838544172e-08*(((330.0-6.0*a)*a-7920.0)*bsq+a*(a*(a*((110.0-1.0*a)*a-7920.0)+332640.0)-6652800.0)+39916800.0),
	  2.08767569878681e-09*b*((132.0-4.0*a)*bsq+a*(a*((528.0-5.0*a)*a-35640.0)+1330560.0)-19958400.0));

  LatticeComplex a2 = 
    cmplx(2.08767569878681e-09*(((396.0-6.0*a)*a-11880.0)*bsq+a*(a*(a*((132.0-1.0*a)*a-11880.0)+665280.0)-19958400.0)+239500800.0),
	  2.505210838544172e-08*b*(bsq+a*(a*(4.0*a-330.0)+15840.0)-332640.0));


  LatticeColorMatrix aux2 = a2*msq + a1*m + a0;
    
  if ( cflag == REUNITARIZE_LABEL )
  {
    int numbad;
    reunit(aux2, numbad, REUNITARIZE_LABEL);
    if ( numbad > 0 )
    {
      QDPIO::cerr << "found " << numbad << " some matrices violating unitarity" << endl;
      QDP_abort(1);
/*#   eesun (a, aux1, aa, bb, aux2); */
    }
  }
  a = aux2;
  
  END_CODE("expsu3");
}
