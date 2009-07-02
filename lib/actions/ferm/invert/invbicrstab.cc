// $Id: invbicrstab.cc,v 3.1 2009-07-02 22:11:03 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invbicrstab.h"

namespace Chroma {

  template<typename T, typename CR>
SystemSolverResults_t
InvBiCRStab_a(const LinearOperator<T>& A,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      int MaxBiCGStab, 
	      enum PlusMinus isign)

{
  SystemSolverResults_t ret;
  StopWatch swatch;
  FlopCounter flopcount;
  flopcount.reset();
  const Subset& s = A.subset();
  bool convP = false;

  swatch.reset();
  swatch.start();

  Double chi_sq =  norm2(chi,s);
  flopcount.addSiteFlops(4*Nc*Ns,s);


  Double rsd_sq =  RsdBiCGStab*RsdBiCGStab*chi_sq;

  // First get r = r0 = chi - A psi
  T r;
  T r0;

  // Get A psi, use r0 as a temporary
  A(r0, psi, isign);
  flopcount.addFlops(A.nFlops());

  // now work out r= chi - Apsi = chi - r0
  r[s] = chi - r0;
  flopcount.addSiteFlops(2*Nc*Ns,s);


  // The main difference between BICGStab and BiCRStab
  // The shadow residual r0* -> A^\dagger r_0*
#if 1
  if( isign == PLUS ) { 
    A(r0,r, MINUS);
  }
  else { 
    A(r0,r, PLUS);
  }
#else
  A(r0,r,isign);
#endif

  // Everything else stays the same
 
  // Now initialise v = p = 0
  T p;
  T v;

  p[s] = zero;
  v[s] = zero;

  T tmp;
  T t;

  ComplexD rho, rho_prev, alpha, omega;

  // rho_0 := alpha := omega = 1
  // Iterations start at k=1, so rho_0 is in rho_prev
  rho_prev = Double(1);
  alpha = Double(1);
  omega = Double(1);

  // The iterations 
  for(int k = 1; k <= MaxBiCGStab && !convP ; k++) { 
    
    // rho_{k+1} = < r_0 | r >
    rho = innerProduct(r0,r,s);


    if( toBool( real(rho) == 0 ) && toBool( imag(rho) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: rho = 0" << endl;
      QDP_abort(1);
    }

    // beta = ( rho_{k+1}/rho_{k})(alpha/omega)
    ComplexD beta;
    beta = ( rho / rho_prev ) * (alpha/omega);
    
    // p = r + beta(p - omega v)

    // first work out p - omega v 
    // into tmp
    // then do p = r + beta tmp
    CR omega_r = omega;
    CR beta_r = beta;
    tmp[s] = p - omega_r*v;
    p[s] = r + beta_r*tmp;


    // v = Ap
    A(v,p,isign);


    // alpha = rho_{k+1} / < r_0 | v >
    // put <r_0 | v > into tmp
    DComplex ctmp = innerProduct(r0,v,s);


    if( toBool( real(ctmp) == 0 ) && toBool( imag(ctmp) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: <r_0|v> = 0" << endl;
      QDP_abort(1);
    }

    alpha = rho / ctmp;

    // Done with rho now, so save it into rho_prev
    rho_prev = rho;

    // s = r - alpha v
    // I can overlap s with r, because I recompute it at the end.
    CR alpha_r = alpha;
    r[s]  -=  alpha_r*v;


    // t = As  = Ar 
    A(t,r,isign);
    // omega = < t | s > / < t | t > = < t | r > / norm2(t);

    // This does the full 5D norm
    Double t_norm = norm2(t,s);


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
    omega_r = omega;
    alpha_r = alpha;
    tmp[s] = psi + omega_r*r;   
    psi[s] = tmp + alpha_r*p;



    // r = s - omega t = r - omega t1G

    
    r[s] -= omega_r*t;


    Double r_norm = norm2(r,s);


    //    QDPIO::cout << "Iteration " << k << " : r = " << r_norm << endl;
    if( toBool(r_norm < rsd_sq ) ) {
      convP = true;
      ret.resid = sqrt(r_norm);
      ret.n_count = k;

    }
    else { 
      convP = false;
    }

    //-------BiCGStab Flopcounting --------------------------------------
    // flopcount.addSiteFlops(8*Nc*Ns,s);     // <r0|r>
    // flopcount.addSiteFlops(16*Nc*Ns,s);    // p = r + beta p - beta_omega v
    // flopcount.addSiteFlops(8*Nc*Ns,s);  //  <r0 | v>
    // flopcount.addSiteFlops(8*Nc*Ns,s);  //  r -= alpha v
    // flopcount.addSiteFlops(8*Nc*Ns, s); //  < t, r>
    // flopcount.addSiteFlops(4*Nc*Ns, s); //  < t, t>     
    // flopcount.addSiteFlops(16*Nc*Ns,s); // psi += omega r + alpha_p
    // flopcount.addSiteFlops(8*Nc*Ns,s); // r -=omega t
    // flopcount.addSiteFlops(4*Nc*Ns,s); // norm2(r)
    // flopcount.addFlops(2*A.nFlops());  // = 80*Nc*Ns cbsite flops + 2*A
    //----------------------------------------------------------------------
    flopcount.addSiteFlops(80*Nc*Ns,s);
    flopcount.addFlops(2*A.nFlops());


  }
  
  swatch.stop();

  QDPIO::cout << "InvBiCRStab: k = " << ret.n_count << " resid = " << ret.resid << endl;
  flopcount.report("invbicrstab", swatch.getTimeInSeconds());

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
  return InvBiCGStab_a<LatticeFermion, Complex>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}
#endif

template<>
SystemSolverResults_t
InvBiCRStab(const LinearOperator<LatticeFermionF>& A,
	    const LatticeFermionF& chi,
	    LatticeFermionF& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCRStab_a<LatticeFermionF, ComplexF>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

template<>
SystemSolverResults_t
InvBiCRStab(const LinearOperator<LatticeFermionD>& A,
	    const LatticeFermionD& chi,
	    LatticeFermionD& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCRStab_a<LatticeFermionD, ComplexD>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

// Staggered
template<>
SystemSolverResults_t
InvBiCRStab(const LinearOperator<LatticeStaggeredFermion>& A,
	    const LatticeStaggeredFermion& chi,
	    LatticeStaggeredFermion& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCRStab_a<LatticeStaggeredFermion, Complex>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

}  // end namespace Chroma
