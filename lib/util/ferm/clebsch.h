// -*- C++ -*-
/*! \file
 * \brief CG routines
 */

#ifndef __clebsch_h__
#define __clebsch_h__

#include <complex>

namespace Hadron
{
  /** @brief  Calculates Clebsch-Gordon coeficents.
   *
   * <b> Returns </b>
   *
   *  \f$\left(j_1,m_1;j_2,m_2|J,M\right) \f$
   *
   * NOTE: the args j1, m1, etc. are 2*spin
   */
  double clebsch(int j1,int m1,int j2,int m2,int J,int M);


  /** Calculates the Wigner d-functions.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   *
   * \returns \f$d^{j}_{m,n}(\beta)\f$.
   *
   */
  double Wigner_d(int J,int M,int N,double __beta);


  /** Calculates the Wigner D-functions.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   *
   * \returns \f$D^{j}_{m,n}(\alpha, \beta, \gamma)\f$.
   *
   */
  std::complex<double> Wigner_D(int J,int M,int N,double __alpha, double __beta, double __gamma);


  /** Calculates sum over product of two Wigner D-functions.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   *
   * \returns \f$\sum_{mm} D^{j}_{m,mm}(\alpha1, \beta1, \gamma1) D^{j}_{mm,n}(\alpha2, \beta2, \gamma2)\f$.
   *
   */
  std::complex<double> Wigner_D_2rot(int J,int M,int N, double __alpha1, double __beta1, double __gamma1,
				     double __alpha2, double __beta2, double __gamma2);


  /** @brief  Calculates Clebsch-Gordon coeficents.
   *
   * <b> Returns </b>
   *
   *  \f$\left(j_1,m_1;j_2,m_2|J,M\right) \f$
   *
   * NOTE: the args j1, m1, etc. and are 2*spin.
   * This routine is 1-based, where m=+j is mapped to row 1, then m=j-1 to row 2, etc.
   */
  double Clebsch1O(int j1,int m1,int j2,int m2,int J,int M);


}  // end namespace Hadron

#endif
