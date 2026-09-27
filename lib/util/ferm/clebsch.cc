/*! \file
 * \brief CG routines
 */

#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "clebsch.h"

using namespace std;

// define a few macros
#define MAX(x,y) (x>y ? x : y)
#define MIN(x,y) (x<y ? x : y)

namespace Hadron
{

  //_____________________________________________________________________________

  /// Utility function used by Clebsch
  inline double dfact(double __x){
    if((__x < 0.00001) && (__x >= 0.0)) return 1.;
    if(__x < 0) return 0.;
    return __x*dfact(__x - 1.);
  }
  //_____________________________________________________________________________

  /// Returns i! - but gives junk for i >12
  /*  inline int factorial(int __i) {
    int f = 1;
    if((__i == 0)||(__i == 1)) f = 1;
    else{
      while(__i > 0){
	f = f*__i;
	__i--;
      }
    }
    return f;
    } */

  // Standard-library factorial; no external analysis-library dependency.
  inline double factorial(int i){
    return i < 0 ? 0.0 : std::tgamma(double(i) + 1.0);
  }


  //_____________________________________________________________________________

  /** @brief  Calcultates Clebsch-Gordon coeficents.
   *
   * <b> Returns </b>
   *
   *  \f$\left(j_1,m_1;j_2,m_2|J,M\right) \f$
   *
   * Note: This function was copied from one written by Denis Weygand.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   */
  double clebsch(int j1,int m1,int j2,int m2,int J,int M)
  {
    if((m1 + m2) != M) return 0.;

    // convert to pure integers (each 2*spin)
//  int j1 = (int)(2.*__j1);
//  int m1 = (int)(2.*__m1);
//  int j2 = (int)(2.*__j2);
//  int m2 = (int)(2.*__m2);
//  int J = (int)(2.*__J);
//  int M = (int)(2.*__M);

    double n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,A,exp;
    int nu = 0;

    double sum = 0;
    while(((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 )) { nu++;}
    while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0)
	   &&((n1=(j2+J+m1)/2-nu) >= 0 )){
      d3=((j1-j2-M)/2+nu);
      n2=((j1-m1)/2+nu);
      d0=dfact((double) nu);
      exp=nu+(j2+m2)/2;
      n0 = (double) pow(-1.,exp);
      sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
      nu++;
    }

    if (sum == 0) return 0;

    n0 = J+1;
    n1 = dfact((double) (J+j1-j2)/2);
    n2 = dfact((double) (J-j1+j2)/2);
    n3 = dfact((double) (j1+j2-J)/2);
    n4 = dfact((double) (J+M)/2);
    n5 = dfact((J-M)/2);

    d0 = dfact((double) (j1+j2+J)/2+1);
    d1 = dfact((double) (j1-m1)/2);
    d2 = dfact((double) (j1+m1)/2);
    d3 = dfact((double) (j2-m2)/2);
    d4 = dfact((double) (j2+m2)/2);

    A = ((double) (n0*n1*n2*n3*n4*n5))/((double) (d0*d1*d2*d3*d4));

    return sqrt(A)*sum;
  }
  //_____________________________________________________________________________

  /** Calculates the Wigner d-functions.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   *
   * \returns \f$d^{j}_{m,n}(\beta)\f$.
   *
   */
  double Wigner_d(int J,int M,int N,double __beta)
  {
//  int J = (int)(2.*__j);
//  int M = (int)(2.*__m);
//  int N = (int)(2.*__n);
    int temp_M, k, k_low, k_hi;
    double const_term = 0.0, sum_term = 0.0, d = 1.0;
    int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
    int kmn1, kmn2, jmnk, jmk, jnk;
    double kk;

    if (J < 0 || abs (M) > J || abs (N) > J) {
      cerr << endl;
      cerr << "d: you have entered an illegal number for J, M, N." << endl;
      cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J."
	   << endl;
      cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
      return 0.;
    }

    if (__beta < 0) {
      __beta = fabs (__beta);
      temp_M = M;
      M = N;
      N = temp_M;
    }

    m_p_n = (M + N) / 2;
    j_p_m = (J + M) / 2;
    j_m_m = (J - M) / 2;
    j_p_n = (J + N) / 2;
    j_m_n = (J - N) / 2;

    kk = (double)factorial(j_p_m)*(double)factorial(j_m_m)
      *(double)factorial(j_p_n) * (double)factorial(j_m_n) ;
    const_term = pow((-1.0),(j_p_m)) * sqrt(kk);

    k_low = MAX(0, m_p_n);
    k_hi = MIN(j_p_m, j_p_n);

    for (k = k_low; k <= k_hi; k++)
    {
      kmn1 = 2 * k - (M + N) / 2;
      jmnk = J + (M + N) / 2 - 2 * k;
      jmk = (J + M) / 2 - k;
      jnk = (J + N) / 2 - k;
      kmn2 = k - (M + N) / 2;

      sum_term += pow ((-1.0), (k)) *
	((pow (cos (__beta / 2.0), kmn1)) * (pow (sin (__beta / 2.0), jmnk))) /
	(factorial (k) * factorial (jmk) * factorial (jnk) * factorial (kmn2));
    }

    d = const_term * sum_term;
    return d;
  }
  //_____________________________________________________________________________



  /** Calculates the Wigner D-functions.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   *
   * \returns \f$D^{j}_{m,n}(\alpha, \beta, \gamma)\f$.
   *
   */
  std::complex<double> Wigner_D(int J,int M,int N,double __alpha, double __beta, double __gamma)
  {
    const std::complex<double> complex_i(0.0, 1.0);

    if (J < 0 || abs (M) > J || abs (N) > J) {
      cerr << endl;
      cerr << "Wigner_D: you have entered an illegal number for J, M, N." << endl;
      cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J."
	   << endl;
      cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
      exit(1);
      return 0.;
    }

    double expPhase = - double(__alpha) * double(M) / double(2.0) - double(__gamma) * double(N) / double(2.0);
    auto D = cos(expPhase) + complex_i * sin(expPhase);
    D = D * Wigner_d(J,M,N,__beta);

    // Debugging
//    std::cout << __func__ << "(" << J << "," << M << "," << N << "; " << __alpha << "," << __beta << "," << __gamma << ") = ";
//    std::cout << D << std::endl;

    return D;
  }
  //_____________________________________________________________________________


  /** Calculates sum over product of two Wigner D-functions.
   *
   * NOTE also: the args j1, m1, etc. are 2*spin,
   *
   * \returns \f$\sum_{mm} D^{j}_{m,mm}(\alpha1, \beta1, \gamma1) D^{j}_{mm,n}(\alpha2, \beta2, \gamma2)\f$.
   *
   */
  std::complex<double> Wigner_D_2rot(int J,int M,int N, double __alpha1, double __beta1, double __gamma1,
				     double __alpha2, double __beta2, double __gamma2)
  {
    std::complex<double> D = 0.0;
    for (int MM = -J; MM <= J; MM=MM+2)
    {
      D += Wigner_D(J,M,MM, __alpha1,__beta1,__gamma1) * Wigner_D(J,MM,N, __alpha2,__beta2,__gamma2);
    }
    return D;
  }
  //_____________________________________________________________________________


  namespace
  {
    int rtM(int twoJ, int r)
    {
      if (r < 1 || r > twoJ+1)
      {
	std::cerr << __func__ <<  ": r out of bounds: 1 <= " << r << " <= " << twoJ+1 << "\n";
	exit(1);
      }

      return -2*(r - 1) + twoJ;
    }
  }

  /** Calculates Clebsch-Gordon coeficents.
   *
   * <b> Returns </b>
   *
   *  \f$\left(j_1,m_1;j_2,m_2|J,M\right) \f$
   *
   * NOTE: the args j1, m1, etc. and are 2*spin.
   * This routine is 1-based, where m=+j is mapped to row 1, then m=j-1 to row 2, etc.
   */
  double Clebsch1O(int j1,int r1,int j2,int r2,int J,int r)
  {
    return clebsch(j1, rtM(j1,r1), j2, rtM(j2,r2), J, rtM(J,r));
  }

} // namespace Hadron
