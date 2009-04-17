// $Id: invbicgstab.cc,v 3.3 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invbicgstab.h"

namespace Chroma {

template<typename T>
SystemSolverResults_t
InvBiCGStab_a(const LinearOperator<T>& A,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      int MaxBiCGStab, 
	      enum PlusMinus isign)

{
  SystemSolverResults_t ret;

  const Subset& s = A.subset();

  bool convP = false;
  Real chi_sq =  Real(norm2(chi,s));
  Real rsd_sq =  RsdBiCGStab*RsdBiCGStab*chi_sq;

  // First get r = r0 = chi - A psi
  T r;
  T r0;

  // Get A psi, use r0 as a temporary
  A(r0, psi, isign);

  // now work out r= chi - Apsi = chi - r0
  r[s] = chi - r0;

  // Also copy back to r0. We are no longer in need of the
  // nth component
  r0[s] = r;
  
  // Now we have r = r0 = chi - Mpsi
 
  // Now initialise v = p = 0
  T p;
  T v;

  p[s] = zero;
  v[s] = zero;

  T tmp;
  T t;

  Complex rho, rho_prev, alpha, omega;

  // rho_0 := alpha := omega = 1
  // Iterations start at k=1, so rho_0 is in rho_prev
  rho_prev = Real(1);
  alpha = Real(1);
  omega = Real(1);

  // The iterations 
  for(int k = 1; k <= MaxBiCGStab && !convP ; k++) { 
    
    // rho_{k+1} = < r_0 | r >
    rho = innerProduct(r0,r,s);

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
    tmp[s] = p - omega*v;
    p[s] = r + beta*tmp;
    
    // v = Ap
    A(v,p,isign);

    // alpha = rho_{k+1} / < r_0 | v >
    // put <r_0 | v > into tmp
    Complex ctmp = innerProduct(r0,v,s);

    if( toBool( real(ctmp) == 0 ) && toBool( imag(ctmp) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: <r_0|v> = 0" << endl;
      QDP_abort(1);
    }

    alpha = rho / ctmp;

    // Done with rho now, so save it into rho_prev
    rho_prev = rho;

    // s = r - alpha v
    // I can overlap s with r, because I recompute it at the end.
    r[s]  -=  alpha*v;
    
    // t = As  = Ar 
    A(t,r,isign);
    
    // omega = < t | s > / < t | t > = < t | r > / norm2(t);

    // This does the full 5D norm
    Real t_norm = norm2(t,s);
    if( toBool(t_norm == 0) ) { 
      QDPIO::cerr << "Breakdown || Ms || = || t || = 0 " << endl;
      QDP_abort(1);
    }

    // accumulate <t | s > = <t | r> into omega
    omega = innerProduct(t,r,s);
    omega /= t_norm;

    // psi = psi + omega s + alpha p 
    //     = psi + omega r + alpha p
    //
    // use tmp to compute psi + omega r
    // then add in the alpha p
    tmp[s] = psi + omega*r;
    psi[s] = tmp + alpha*p;


    // r = s - omega t = r - omega t1G
    Double r_norm;
    
    r[s] -= omega*t;

    r_norm = norm2(r,s);

    //    QDPIO::cout << "Iteration " << k << " : r = " << r_norm << endl;
    if( toBool(r_norm < rsd_sq ) ) {
      convP = true;
      ret.resid = sqrt(r_norm);
      ret.n_count = k;

    }
    else { 
      convP = false;
    }



  }

  QDPIO::cout << "InvBiCGStab: k = " << ret.n_count << " resid = " << ret.resid << endl;
  if ( ret.n_count == MaxBiCGStab ) { 
    QDPIO::cerr << "Nonconvergence of BiCGStab. MaxIters reached " << endl;
  }

  return ret;
}

#if 0
// Fix here for now
template<>
SystemSolverResults_t
InvBiCGStab(const LinearOperator<LatticeFermion>& A,
	    const LatticeFermion& chi,
	    LatticeFermion& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCGStab_a(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}
#endif

template<>
SystemSolverResults_t
InvBiCGStab(const LinearOperator<LatticeFermionF>& A,
	    const LatticeFermionF& chi,
	    LatticeFermionF& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCGStab_a(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

template<>
SystemSolverResults_t
InvBiCGStab(const LinearOperator<LatticeFermionD>& A,
	    const LatticeFermionD& chi,
	    LatticeFermionD& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCGStab_a(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

// Staggered
template<>
SystemSolverResults_t
InvBiCGStab(const LinearOperator<LatticeStaggeredFermion>& A,
	    const LatticeStaggeredFermion& chi,
	    LatticeStaggeredFermion& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCGStab_a(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

}  // end namespace Chroma
