// $Id: eesu3.cc,v 1.3 2004-03-03 01:50:40 edwards Exp $
/*! \file
 *  \brief Exactly exponentiate a SU(3) lie algebra element
 */

#include "chromabase.h"
#include "util/gauge/eesu3.h"

using namespace QDP;

//! Exactly exponentiate a SU(3) lie algebra element
/*!
 * \ingroup gauge
 *
 * Follows hep-lat/0311018 . Convention change, there C&M expect a Q
 * which is traceless HERMITIAN. Here, we expect traceless ANTI-HERMITIAN.
 *
 *  \param m       Initially a traceless/anti-herm. matrix      (Modify)
 *
 */

void eesun(LatticeColorMatrix3& m)
{
  START_CODE("eesu3");

  // Goal, evaluate
  // exp(i*Q) = f0*I + f1*Q + f2*Q^2

  // NOTE, pass m = i*Q  ->  Q = -i*m

  // Q = -i*m
  LatticeColorMatrix3 Q  = timesMinusI(m);
  
  // Q2 = Q*Q
  LatticeColorMatrix3 Q2 = Q*Q;
  
  // c0 = tr(Q^3)/3 = real(tr(Q^3))/3
  Real third = 1.0/3.0;
  LatticeReal c0 = third*real(trace(Q*Q2));
  
  // c1 = tr(Q*Q)/2
  LatticeReal c1 = Real(0.5)*real(trace(Q2));
  
  // c0max = 2 * (c1/3)^{3/2}
//exp(Real(1.5)*log(third*c1));
  LatticeReal c0max = 2*pow(third*c1, Real(1.5));

  // Basic terms
  // NOTE: take ABSOLUTE value of c0 and later use relation
  // that  f_j(-c0,c1) = (-1)^{j}*conj(f_j(c0,c1))
  LatticeReal theta = acos(abs(c0)/c0max);
  LatticeReal u = sqrt(third*c1)*cos(third*theta);
  LatticeReal w = sqrt(c1)*sin(third*theta);

  // Carefully control  xi0
  LatticeReal w2  = w*w;
  LatticeReal xi0 = where(abs(w) < 0.05, 
			  1 - Real(1/6)*w2*(1-Real(1/20)*w2*(1-Real(1/42)*w2)),
			  sin(w)/w);

  // The h_i
  LatticeComplex e2u = cmplx(cos(2*u),sin(2*u));
  LatticeComplex emu = cmplx(cos(u),-sin(u));
  LatticeReal u2     = u*u;
  LatticeReal cosw   = cos(w);

  LatticeComplex h0 = (u2-w2)*e2u + emu*(8*u2*cosw + timesI(2*u*(3*u2+w2)*xi0));
  LatticeComplex h1 = 2*u*e2u - emu*(2*u*cosw - timesI((3*u2-w2)*xi0));
  LatticeComplex h2 = e2u - emu*(cosw + timesI(3*u*xi0));

  //
  // The f_j . Note, by using abs(c0) above, we do not have
  // a singularity problem. So, correct for abs(c0) by using
  // f_j(-c0,c1) = (-1)^{j}*conj(f_j(c0,c1))
  //
  LatticeReal    fact = Real(1.0) / (9*u2 - w2);
  LatticeComplex f0   = where(c0 > 0, h0*fact,  conj(h0*fact));
  LatticeComplex f1   = where(c0 > 0, h1*fact, -conj(h1*fact));
  LatticeComplex f2   = where(c0 > 0, h2*fact,  conj(h2*fact));

  // Finally
  m = f0 + f1*Q + f2*Q2;

  END_CODE("eesu3");
}
