// $Id: invbicgstab.cc,v 3.5 2009-06-30 15:52:10 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invbicgstab.h"


#include "actions/ferm/invert/bicgstab_kernels.h"

using namespace Chroma::BiCGStabKernels;

namespace Chroma {
  template<typename T, typename CR>
SystemSolverResults_t
InvBiCGStab_a(const LinearOperator<T>& A,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      int MaxBiCGStab, 
	      enum PlusMinus isign)

{
  
  SystemSolverResults_t ret;
  StopWatch swatch;
  swatch.reset();
  swatch.start();

  initKernels();
  
  FlopCounter flopcount;
  flopcount.reset();
  const Subset& sub = A.subset();
  bool convP = false;


  T r0,r,u,f0,v,z,q,t,s;
  DComplex sigma_n_2, sigma_n_1, sigma_n, pi_n_1, pi_n, phi_n_1,phi_n;
  DComplex tau_n_1, tau_n, rho_n_1, rho_n, alpha_n_1, alpha_n, omega_n_1, omega_n;
  DComplex theta_n, eta_n, gamma_n, beta_n, delta_n;
  Double   kappa_n,rnorm_prev;

  DComplex tcmpx1, tcmpx2; // Temporaries
  T tvec1;

  Double chi_sq =  norm2(chi,sub);
  flopcount.addSiteFlops(4*Nc*Ns,sub);
  Double rsd_sq =  RsdBiCGStab*RsdBiCGStab*chi_sq;




  // First get r = r0 = chi - A psi
  // Get A psi, use r0 as a temporary
  A(r0, psi, isign);
  flopcount.addFlops(A.nFlops());

  // now work out r= chi - Apsi = chi - r0
  r[sub] = chi - r0;
  flopcount.addSiteFlops(2*Nc*Ns,sub);


  // Also copy back to r0. We are no longer in need of the
  // nth component
  r0[sub] = r;
  
  // u = Ar;
  A(u, r, isign);
  flopcount.addFlops(A.nFlops());


  // f0 = A^\dag r
  if( isign == PLUS ) { 
    A(f0, r, MINUS);
  }
  else { 
    A(f0, r, PLUS);
  }
  flopcount.addFlops(A.nFlops());

  // q_n = v_n = z_n = 0;
  q[sub]=zero;
  v[sub]=zero;
  z[sub]=zero;

  sigma_n_2 = Double(0);
  pi_n_1 = Double(0);
  tau_n_1 = Double(0);

  // Deviate from paper. Paper says phi_n = 0. Follow code from PETSc 
  phi_n_1 = innerProduct(r0,r0,sub);
  sigma_n_1 = innerProduct(r0,u, sub);
  flopcount.addSiteFlops(16*Nc*Ns,sub);   // each inner product is 8

  rho_n_1 = Double(1);
  alpha_n_1 = Double(1);
  omega_n_1 = Double(1);

  // Main Loop starts here
  for(int k = 1; k <= MaxBiCGStab && !convP ; k++) { 

    // 16 Flops: NB I pulled the -omega factor out. This may
    // or may not be good.
    rho_n = phi_n_1 - omega_n_1*(sigma_n_2 - alpha_n_1*pi_n_1);
    

    if( k == 1)  {
      delta_n = rho_n;
    }
    else { 
      /* NB:  The paper says rho_n/rho_n_1, rather than rho_n/tau_n_1 
       * Again here for now I will defer to the petsc code. */
      delta_n = rho_n /tau_n_1;  // 11 Flops
    }

    beta_n = delta_n /omega_n_1; // 11 Flops

    tau_n = sigma_n_1 + beta_n*tau_n_1 - delta_n*pi_n_1; // 16 flops


    alpha_n = rho_n/tau_n;       // 11 Flops


    /* Now some vector update.
     * 
     * PETSc claims code in paper is wrong for z_n update
     *
     * z_n = alpha_n*r_{n-1} + (alpha_n/alpha_{n-1})*beta_n z_{n-1}
     *         -alpha_n*delta_n*v_{n-1}
     *
     * v_n = u_{n-1] + beta_n*v_{n-1} - delta_n*q_n{-1}
     */

    // z_n 
    tcmpx1 = (alpha_n/alpha_n_1)*beta_n; // 11 Flops for the ratio, 6 for mul.
    tcmpx2 = alpha_n*delta_n;            // 6 Flops

    // ZV Updates
    ibicgstab_zvupdates(r,z,v,u,q, alpha_n, tcmpx1,tcmpx2,beta_n,delta_n, sub);

#if 0
    tvec1[sub] = tcmpx1*z;
    z[sub] = alpha_n*r + tvec1;
    z[sub] -= tcmpx2*v;
#endif



#if 0
    tvec1[sub] = v;
    v[sub] = u + beta_n*tvec1;
    v[sub]-= delta_n*q;
#endif


    

    /* q = A*v */
    A(q,v,isign);

    /* Does all of:
     *
     *  s = r-alpha_n*v_n
     *  t = u - alpha_n q_n 
     *  phi_n = < r0, s>
     *  gamma_n = <f0, s>
     *  pi      = <r0, q>
     *  eta     = <f0, t>
     *  theta   = <s,  t>
     *  kappa   = || t ||^2
     *  rnorm_prev = || r ||^2 
     * 
     * in principle with 1 big allreduce */
    ibicgstab_stupdates_reduces(alpha_n,
				r, u, v, q, r0, f0,
				s, t, 
				phi_n,pi_n,gamma_n,eta_n,theta_n,kappa_n,rnorm_prev, sub);
    
#if 0
    s[sub] = r - alpha_n*v;
    t[sub] = u - alpha_n*q;

    /* Now 5 inner products + 1 norm + the residual norm from last iteration */
    phi_n = innerProduct(r0,s,sub);
    gamma_n = innerProduct(f0,s,sub);

    pi_n  = innerProduct(r0,q,sub);
    eta_n   = innerProduct(f0,t,sub);
    theta_n = innerProduct(s,t,sub);
    kappa_n  = norm2(t,sub);
    rnorm_prev = norm2(r,sub);
#endif

    // 5 inner products with 8Fl each=> 40, 2 norms with 4 Fl each = 8
    flopcount.addSiteFlops( 54*Nc*Ns, sub);
    flopcount.addFlops(A.nFlops()+88); // not first iteration delta more comlex
    
    // Check whether previous iteration converged:
    if( toBool( rnorm_prev <= rsd_sq ) ) { 
      convP = true;
      ret.n_count = k; // Should this be k-1?
      ret.resid = sqrt(rnorm_prev);
    }
    else { 
      // Not converged: 
      omega_n = theta_n/ kappa_n;           // 11 Flops
      sigma_n = gamma_n - omega_n*eta_n;    // 8 Flops

      // r = s - omega t
      // x = x + z + omega s
      ibicgstab_rxupdate(omega_n,
			    s,
			    t,
			    z,
			    r, 
			    psi,
			    sub);
				

#if 0
      // Update r
      r[sub] = s - omega_n*t;               // 8 Flops/complex

      // Update x 
      psi[sub] += z;                        // 2 flops/complex
      psi[sub] += omega_n*s;                // 8 flops/complex
#endif

      // u = A r
      A(u,r, isign);


      flopcount.addSiteFlops( 18*Nc*Ns, sub);      
      flopcount.addFlops(A.nFlops()+19);

      // Update n_1 locations with n locations
      sigma_n_2 = sigma_n_1;
      sigma_n_1 = sigma_n;
      pi_n_1 = pi_n;
      phi_n_1 = phi_n;
      alpha_n_1 = alpha_n;
      tau_n_1 = tau_n;
      rho_n_1 = rho_n;
      omega_n_1 = omega_n;
    }
    

  }



  swatch.stop();
  if( ret.n_count > 1 ) flopcount.addFlops(-11); // For first iter

  QDPIO::cout << "InvBiCGStab: k = " << ret.n_count << " resid = " << ret.resid << endl;
  flopcount.report("invbicgstab", swatch.getTimeInSeconds());

  if ( ret.n_count == MaxBiCGStab ) { 
    QDPIO::cerr << "Nonconvergence of BiCGStab. MaxIters reached " << endl;
  }

  finishKernels();
  return ret;
}


#if 0
  template<typename T, typename CR>
SystemSolverResults_t
InvBiCGStab_a(const LinearOperator<T>& A,
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

  QDPIO::cout << "InvBiCGStab: k = " << ret.n_count << " resid = " << ret.resid << endl;
  flopcount.report("invbicgstab", swatch.getTimeInSeconds());

  if ( ret.n_count == MaxBiCGStab ) { 
    QDPIO::cerr << "Nonconvergence of BiCGStab. MaxIters reached " << endl;
  }

  return ret;
}
#endif
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
InvBiCGStab(const LinearOperator<LatticeFermionF>& A,
	    const LatticeFermionF& chi,
	    LatticeFermionF& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvBiCGStab_a<LatticeFermionF, ComplexF>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
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
  return InvBiCGStab_a<LatticeFermionD, ComplexD>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
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
  return InvBiCGStab_a<LatticeStaggeredFermion, Complex>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

}  // end namespace Chroma
