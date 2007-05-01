// $Id: invbicgstab_array.cc,v 3.3 2007-05-01 12:50:12 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invbicgstab_array.h"

namespace Chroma {

template<typename T>
SystemSolverResults_t
InvBiCGStab_a(const LinearOperatorArray<T>& A,
	      const multi1d<T> & chi,
	      multi1d<T>& psi,
	      const Real& RsdBiCGStab, 
	      int MaxBiCGStab) 
{
  SystemSolverResults_t ret;

  const int N = psi.size();

  const Subset& s = A.subset();

  // Not converged
  bool convP = false;

  // || chi ||^2
  Real chi_sq = norm2(chi,s);

  // Target residuum ( epsilon*|| chi || )^2
  Real rsd_sq =  RsdBiCGStab*RsdBiCGStab*chi_sq;

  // First get r = r0 = chi - A psi
  multi1d<T> r(N);
  multi1d<T> r0(N);

  // Get A psi, use r0 as a temporary
  A(r0, psi, PLUS);

  // now work out r= chi - Apsi = chi - r0
  for(int n=0; n < N; n++) { 
    r[n][s] = chi[n] - r0[n];

    // Also copy back to r0. We are no longer in need of the
    // nth component
    r0[n][s] = r[n];
  }

  // Now we have r = r0 = chi - Mpsi
 
  // Now initialise v = p = 0
  multi1d<T> p(N);
  multi1d<T> v(N);

  for(int n = 0; n < N; n++) { 
    p[n][s] = zero;
    v[n][s] = zero;
  }

  // Create a temporary vector for later use
  T tmp;

  // Search vector for t = Ms step.
  multi1d<T> t(N);

  Complex rho, rho_prev, alpha, omega;

  // rho_0 := alpha := omega = 1
  // Iterations start at k=1, so rho_0 is in rho_prev
  rho_prev = Real(1);
  alpha = Real(1);
  omega = Real(1);

  // The iterations 
  for(int k = 1; k <= MaxBiCGStab && !convP ; k++) { 
    
    // rho_{k+1} = < r_0 | r >
    rho = Double(0);
    for(int n=0; n < N; n++) {
      rho += innerProduct(r0[n],r[n],s);
    }

    
    if( toBool( real(rho) == 0 ) && toBool( imag(rho) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: rho = 0" << endl;
      QDP_abort(1);
    }

    // beta = ( rho_{k+1}/rho_{k})(alpha/omega)
    Complex beta;
    beta = ( rho / rho_prev ) * (alpha/omega);
    
    
    // p = r + beta(p - omega v)

    // first work out p - omega v 
    // into tmp
    // then do p = r + beta tmp
    for(int n=0; n < N; n++) {
      tmp[s] = p[n] - omega*v[n];
      p[n][s] = r[n] + beta*tmp;
    }

    // v = Ap
    A(v,p,PLUS);

    // alpha = rho_{k+1} / < r_0 | v >

    // put <r_0 | v > into ctmp
    Complex ctmp=Real(0);
    for(int n = 0; n < N; n++) { 
      ctmp += innerProduct(r0[n],v[n],s);
    }

    if( toBool( real(ctmp) == 0 ) && toBool( imag(ctmp) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: <r_0|v> = 0" << endl;
      QDP_abort(1);
    }

    alpha = rho / ctmp;

    // Done with rho now, so save it into rho_prev
    rho_prev = rho;

    // s = r - alpha v
    // I can overlap s with r, because I recompute it at the end.
    for(int n = 0; n < N; n++) { 
      r[n][s]  -= alpha*v[n];
    }

    // t = As  = Ar 
    A(t,r,PLUS);
    
    // omega = < t | s > / < t | t > = < t | r > / norm2(t);

    // This does the full 5D norm
    Real t_norm = norm2(t,s);

    if( toBool(t_norm == 0) ) { 
      QDPIO::cerr << "Breakdown || Ms || = || t || = 0 " << endl;
      QDP_abort(1);
    }

    // accumulate <t | s > = <t | r> into omega
    omega = Real(0);
    for(int n = 0; n < N; n++) { 
      omega += innerProduct(t[n],r[n],s);
    }
    
    omega /= t_norm;

    // psi = psi + omega s + alpha p 
    //     = psi + omega r + alpha p
    //
    // use tmp to compute psi + omega r
    // then add in the alpha p
    for(int n=0; n < N; n++) {
      tmp[s]    = psi[n] + omega*r[n];
      psi[n][s] = tmp + alpha*p[n];
    }

    // r = s - omega t = r - omega t => r -= omega_t
    for(int n=0; n < N; n++){ 
      r[n][s] -= omega*t[n];
    }
    
    // Check convergence
    Double r_norm = norm2(r, s);

    QDPIO::cout << "Iteration " << k << " : r = " << r_norm << endl;
    if( toBool(r_norm < rsd_sq ) ) {
      convP = true;
      ret.resid = sqrt(r_norm);
      ret.n_count = k;
    }
    else { 
      convP = false;
    }
  }

  if ( ret.n_count == MaxBiCGStab ) { 
    QDPIO::cerr << "Nonconvergence of BiCGStab. MaxIters reached " << endl;
  }

  return ret;
}


// Fix here for now
template<>
SystemSolverResults_t 
InvBiCGStab(const LinearOperatorArray<LatticeFermion>& A,
	    const multi1d<LatticeFermion>& chi,
	    multi1d<LatticeFermion>& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab)


{
  return InvBiCGStab_a(A, chi, psi, RsdBiCGStab, MaxBiCGStab);
}

}  // end namespace Chroma
