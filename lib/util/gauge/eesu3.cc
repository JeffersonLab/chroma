// $Id: eesu3.cc,v 1.1 2004-01-22 22:52:54 edwards Exp $
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
 *  \param m        LatticeColorMatrix          (Modify)
 */

void eesu3(LatticeColorMatrix& m)
{
  LatticeColorMatrix m_sq;
  LatticeColorMatrix aux1;
  LatticeColorMatrix aux2;
  LatticeComplex ap;
  LatticeComplex bp;
  LatticeComplex cp;

  LatticeBoolean singular_delta;
  LatticeBoolean nonsingular_a;
  LatticeBoolean ltmp1;
  LatticeBoolean ltmp2;

  LatticeReal a;
  LatticeReal b;

  LatticeReal c;
  LatticeReal s;
  LatticeReal u;
  LatticeReal gamma;
  LatticeReal phi;
  LatticeReal tmp1;
  LatticeReal tmp2;
  LatticeReal tmp3;
  LatticeReal tmp4;
  LatticeReal alpha;
  LatticeReal beta;
  LatticeReal alphabeta;
  LatticeReal delta1;
  LatticeReal delta2;
  LatticeReal delta3;
  LatticeReal delta;
  LatticeReal invdelta;
  LatticeReal re_e1;
  LatticeReal re_e2;
  LatticeReal re_e3;
  LatticeReal im_e1;
  LatticeReal im_e2;
  LatticeReal im_e3;
  LatticeReal re_ap;
  LatticeReal re_bp;
  LatticeReal re_cp;
  LatticeReal im_ap;
  LatticeReal im_bp;
  LatticeReal im_cp;
  LatticeReal re_ap_tmp;
  LatticeReal re_bp_tmp;
  LatticeReal re_cp_tmp;
  LatticeReal im_ap_tmp;
  LatticeReal im_bp_tmp;
  LatticeReal im_cp_tmp;
  Real dummy;

/*# debug variables */
/*#int num_small_a; */
/*#int num_singular_delta; */
/*#int num_singular_delta1; */
/*#int num_nonsingular; */
/*#int temp; */

  START_CODE("eesu3");

  if ( Nc != 3 )
    QDP_error_exit("can only handle SU[3]");



/*# aux1 = m*m */
  m_sq = m * m;

/*# a = -tr(m*m)/2 */
  a = real(trace(m_sq));
  dummy = 0.5;
  dummy = -dummy;
  a = a * dummy;

/*# aux2 = m*m*m */
  aux1 = m_sq * m;

/*# b = i*tr(m*m*m)/3 = -Im( tr(m*m*m)/3 ) */
  b = imag(trace(aux1));
  dummy = -TO_REAL(1)/TO_REAL(3);
  b = b * dummy;



/*# u = sqrt(a/3); */
  dummy = TO_REAL(1)/TO_REAL(3);
  tmp1 = a * dummy;
  u = sqrt(tmp1);

/*# gamma = b/(2*u**3) */
  dummy = 2;
  tmp1 = u * u;
  tmp2 = tmp1 * u;
  tmp1 = tmp2 * dummy;

  gamma = b / tmp1;

/*# phi = acos(gamma)/3.0E()0 */
  dummy = TO_REAL(1)/TO_REAL(3);
  tmp1 = acos(gamma);
  tmp1 = tmp1 * dummy;

/*# c = cos(phi) */
  c = cos(tmp1);

/*# s = sin(phi) */
  s = sin(tmp1);

/*# alpha = u*(-cos(phi) + sqrt(3)*sin(phi)) */
  tmp1 = -(c * u);
  dummy = sqrt(TO_REAL(3));
  tmp2 = s * dummy;
  tmp2 = tmp2 * u;

  alpha = tmp1;
  alpha += tmp2;

/*# beta = u*(-cos(phi) - sqrt(3)*sin(phi)) */
  beta = tmp1;
  beta -= tmp2;

/*# alphabeta = -2*u*cos(phi) */
  dummy = 2;
  alphabeta = tmp1 * dummy;





/*# delta1 = alpha + 2*beta */
  dummy = 2;
  delta1 = beta * dummy;
  delta1 += alpha;

/*# delta2 = -(2*alpha + beta) */
  delta2 = -(alpha * dummy);
  delta2 -= beta;

/*# delta3 = alpha - beta */
  delta3 = alpha;
  delta3 -= beta;

/*# delta = delta1*delta2*delta3 */
  tmp1 = delta1 * delta2;
  delta = tmp1 * delta3;

/*### debug #### */
/*#temp = vol; */
/*#vol = 4; */
/*#push(nml_out,"Eesu3_info");
  Write(nml_out, a);
  Write(nml_out, b);
  Write(nml_out, alpha);
  Write(nml_out, beta);
  Write(nml_out, delta1);
  Write(nml_out, delta2);
  Write(nml_out, delta3);
  Write(nml_out, delta);
  pop(nml_out); */
/*#vol = temp; */
/*############## */

/*# Start of conditional calculations. */
/*# Check for small tr(m^2). */
  dummy = 1.0E()-04;
  ltmp1 = a < dummy;

/*#num_small_a = CM_global_count_context(ltmp1); */

/*# Make delta one. */

  tmp3 = 1;
  copymask(delta, ltmp1, tmp3, REPLACE);

  tmp1 = alpha * alpha;
  tmp1 += alpha * beta;
  tmp1 += beta * beta;
  tmp2 = alphabeta * alpha;
  tmp2 = tmp2 * beta;
  tmp3 = tmp1 * tmp1;
  tmp4 = tmp2 * tmp1;

/*# ap = 0.5 - (alpha^2+alpha*beta+beta^2)/24 + i*alpha*beta*(alpha+beta)/120 */
  re_ap_tmp = 0.5;
  dummy = TO_REAL(1) / TO_REAL(24);
  re_ap_tmp -= tmp1 * dummy;
  dummy = TO_REAL(1) / TO_REAL(720);
  re_ap_tmp += tmp3 * dummy;
  dummy = TO_REAL(1) / TO_REAL(120);
  im_ap_tmp = tmp2 * dummy;
  dummy = TO_REAL(1) / TO_REAL(2520);
  im_ap -= tmp4 * dummy;
  copymask(re_ap, ltmp1, re_ap_tmp, REPLACE);
  copymask(im_ap, ltmp1, im_ap_tmp, REPLACE);

/*# bp = 1.0 - (alpha^2+alpha*beta+beta^2)/6 + i*alpha*beta*(alpha+beta)/24 */
  re_bp_tmp = 1;
  dummy = TO_REAL(1) / TO_REAL(6);
  re_bp_tmp -= tmp1 * dummy;
  dummy = TO_REAL(1) / TO_REAL(120);
  re_bp_tmp += tmp3 * dummy;
  dummy = TO_REAL(1) / TO_REAL(24);
  im_bp_tmp = tmp2 * dummy;
  dummy = TO_REAL(1) / TO_REAL(360);
  im_bp_tmp -= tmp4 * dummy;
  copymask(re_bp, ltmp1, re_bp_tmp, REPLACE);
  copymask(im_bp, ltmp1, im_bp_tmp, REPLACE);

/*# cp = 1.0 + i*alpha*beta*(alpha+beta)/6 */
  re_cp_tmp = 1;
  dummy = TO_REAL(1) / TO_REAL(6);
  im_cp_tmp = tmp2 * dummy;
  dummy = TO_REAL(1) / TO_REAL(120);
  im_cp_tmp = tmp4 * dummy;
  copymask(re_cp, ltmp1, re_cp_tmp, REPLACE);
  copymask(im_cp, ltmp1, im_cp_tmp, REPLACE);



/*# Check for singularity in delta (the determinant). */
  nonsingular_a = !ltmp1;

  tmp1 = fabs(delta);
  FILL(dummy,FUZZ);
  singular_delta = tmp1 < dummy;
  singular_delta = singular_delta & nonsingular_a;

/*#num_singular_delta = CM_global_count_context(); */

/*# Context is now a singular delta. */
/*# Check if delta1 is singular. If so, set alpha <- -alpha/2 = beta */
  tmp1 = fabs(delta1);
  tmp2 = fabs(delta2);
  tmp3 = fabs(delta3);

  ltmp1 = tmp1 < tmp2;
  ltmp2 = tmp2 < tmp3;
  ltmp1 = ltmp2 & ltmp1;
  ltmp1 = singular_delta & ltmp1;
  copymask(alpha, ltmp1, beta, REPLACE);


/*#num_singular_delta1 = CM_global_count_context(ltmp1); */

/*# Assuming only two degenerate eigenvalues. */
/*# Context is set back to singular_delta */

/*# e1 = exp(i*alpha) */
  re_e1 = cos(alpha);
  im_e1 = sin(alpha);

/*# e2 = exp(-2*i*alpha) */
  dummy = 2;
  dummy = -dummy;
  tmp1 = alpha * dummy;
  re_e2 = cos(tmp1);
  im_e2 = sin(tmp1);

/*# Now, delta <- 9*alpha^2 */
  dummy = TO_REAL(3);
  tmp1 = alpha * dummy;
  delta = tmp1 * tmp1;

/*# ap = ((1-3*i*alpha)*e1 - e2)/delta */
  dummy = TO_REAL(3);
  tmp3 = alpha * dummy;
  re_ap_tmp = re_e1;
  re_ap_tmp += tmp3 * im_e1;
  im_ap_tmp = im_e1;
  im_ap_tmp -= tmp3 * re_e1;

  re_ap_tmp -= re_e2;
  im_ap_tmp -= im_e2;
  copymask(re_ap, singular_delta, re_ap_tmp, REPLACE);
  copymask(im_ap, singular_delta, im_ap_tmp, REPLACE);

/*# bp = ((3*alpha^2-2*i*alpha)*e1 + 2*i*alpha*e2)/delta */
  dummy = TO_REAL(3);
  tmp3 = alpha * dummy;
  tmp3 = tmp3 * alpha;
  re_bp_tmp = tmp3 * re_e1;
  im_bp_tmp = tmp3 * im_e1;
  dummy = 2;
  tmp3 = alpha * dummy;
  re_bp_tmp += tmp3 * im_e1;
  im_bp_tmp -= tmp3 * re_e1;

  re_bp_tmp -= tmp3 * im_e2;
  im_bp_tmp += tmp3 * re_e2;
  copymask(re_bp, singular_delta, re_bp_tmp, REPLACE);
  copymask(im_bp, singular_delta, im_bp_tmp, REPLACE);

/*# cp = alpha^2*((8-6*i*alpha)*e1 + e2)/delta */
  dummy = TO_REAL(8);
  re_cp_tmp = re_e1 * dummy;
  im_cp_tmp = im_e1 * dummy;

  dummy = TO_REAL(6);
  tmp3 = alpha * dummy;
  re_cp_tmp += tmp3 * im_e1;
  im_cp_tmp -= tmp3 * re_e1;

  re_cp_tmp += re_e2;
  im_cp_tmp += im_e2;

  tmp3 = alpha * alpha;
  re_cp_tmp = re_cp_tmp * tmp3;
  im_cp_tmp = im_cp_tmp * tmp3;
  copymask(re_cp, singular_delta, re_cp_tmp, REPLACE);
  copymask(im_cp, singular_delta, im_cp_tmp, REPLACE);



/*# Flip the context and do the case of delta not singular and tr(m^2) */
/*# not small. */
  ltmp1 = !singular_delta;
  ltmp1 = nonsingular_a & ltmp1;

/*#num_nonsingular = CM_global_count_context(); */


/*# e1 = exp(i*alpha) */
  re_e1 = cos(alpha);
  im_e1 = sin(alpha);

/*# e2 = exp(i*beta) */
  re_e2 = cos(beta);
  im_e2 = sin(beta);

/*# e3 = exp(-i*alphabeta) */
  re_e3 = cos(alphabeta);
  im_e3 = sin(alphabeta);

/*# ap = (delta1*e1 + delta2*e2 + delta3*e3)/delta */
  re_ap_tmp = delta1 * re_e1;
  re_ap_tmp += delta2 * re_e2;
  re_ap_tmp += delta3 * re_e3;
  im_ap_tmp = delta1 * im_e1;
  im_ap_tmp += delta2 * im_e2;
  im_ap_tmp += delta3 * im_e3;
  copymask(re_ap, ltmp1, re_ap_tmp, REPLACE);
  copymask(im_ap, ltmp1, im_ap_tmp, REPLACE);

/*# bp = (i*alpha*delta1*e1 + i*beta*delta2*e2 - i*alphabeta*delta3*e3)/delta */
  tmp3 = alpha * delta1;
  re_bp_tmp = -(tmp3 * im_e1);
  im_bp_tmp = tmp3 * re_e1;

  tmp3 = beta * delta2;
  re_bp_tmp -= tmp3 * im_e2;
  im_bp_tmp += tmp3 * re_e2;

  tmp3 = alphabeta * delta3;
  re_bp_tmp += tmp3 * im_e3;
  im_bp_tmp -= tmp3 * re_e3;
  copymask(re_bp, ltmp1, re_bp_tmp, REPLACE);
  copymask(im_bp, ltmp1, im_bp_tmp, REPLACE);

/*# cp = (beta*alphabeta*delta1*e1 + alpha*alphabeta*delta2*e2 */
/*#    - alpha*beta*delta3*e3)/delta */
  tmp3 = beta * alphabeta;
  tmp3 = tmp3 * delta1;
  re_cp_tmp = tmp3 * re_e1;
  im_cp_tmp = tmp3 * im_e1;

  tmp3 = alpha * alphabeta;
  tmp3 = tmp3 * delta2;
  re_cp_tmp += tmp3 * re_e2;
  im_cp_tmp += tmp3 * im_e2;

  tmp3 = alpha * beta;
  tmp3 = tmp3 * delta3;
  re_cp_tmp -= tmp3 * re_e3;
  im_cp_tmp -= tmp3 * im_e3;
  copymask(re_cp, ltmp1, re_cp_tmp, REPLACE);
  copymask(im_cp, ltmp1, im_cp_tmp, REPLACE);


/*# Through with conditional calculations. Restore context and continue. */
/*#push(nml_out,"Context_sizes");
  Write(nml_out, num_small_a);
  Write(nml_out, num_singular_delta);
  Write(nml_out, num_singular_delta1);
  Write(nml_out, num_nonsingular);
  pop(nml_out); */


/*# Calculate 1/delta where delta has been conditionally computed */
  tmp1 = 1;
  invdelta = tmp1 / delta;

  re_ap = re_ap * invdelta;
  im_ap = im_ap * invdelta;
  re_bp = re_bp * invdelta;
  im_bp = im_bp * invdelta;
  re_cp = re_cp * invdelta;
  im_cp = im_cp * invdelta;


/*# Build lattice complex variables out of two lattice real variables. */
  ap = cmplx(re_ap,im_ap);
  bp = cmplx(re_bp,im_bp);
  cp = cmplx(re_cp,im_cp);



/*# Now calculate the exponential */
/*#  exp(M) = ap*M^2 + bp*M + cp*1 */

  aux2 = 1;
  aux1 = aux2 * cp;
  aux1 += m * bp;
  aux1 += m_sq * ap;
  m = aux1;


  END_CODE("eesu3");
}
