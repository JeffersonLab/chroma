// $Id: invibicgstab.cc,v 3.2 2009-09-01 19:51:28 jbulava Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invibicgstab.h"
#include "actions/ferm/invert/bicgstab_kernels.h"

namespace Chroma {

  using namespace BiCGStabKernels;

  template<typename T, typename CR>
SystemSolverResults_t
InvIBiCGStab_a(const LinearOperator<T>& A,
	       const T& chi,
	       T& psi,
	       const Real& RsdBiCGStab,
	       int MaxBiCGStab, 
	       enum PlusMinus isign)

{
  SystemSolverResults_t ret;
  initKernels();
  StopWatch swatch;
  FlopCounter flopcount;
  flopcount.reset();
  const Subset& sub = A.subset();
  bool convP = false;
	ret.n_count = MaxBiCGStab;

  swatch.reset();
  swatch.start();

  Double chi_sq =  norm2(chi,sub);
  flopcount.addSiteFlops(4*Nc*Ns,sub);


  Double rsd_sq =  RsdBiCGStab*RsdBiCGStab*chi_sq;

  // First get r = r0 = chi - A psi
  T r;
  T r0;
  T v;
  T tmp;
  T t;
  T s;
  T u;
  T q;
  T f0; 
  T z;
  ComplexD rhon, rhon_1;             // rho_n, rho_{n-1}
  ComplexD alphan, alphan_1;         // alpha_{n}, alpha_{n-1}
  ComplexD omega;                    // omega_{n} AND omega_{n-1} (before omega update)
  ComplexD taun,taun_1;              // tau_{n}, tau_{n-1}
  RealD kappa;                       
  RealD rnorm_prev;                  // lagged residuum
  ComplexD theta;
  ComplexD sigman_2, sigman_1;       
  ComplexD phi;                      // phi_{n} AND phi_{n-1} before phi update
  ComplexD pi_c;                     // pi_{n} and pi_{n-1} before pi update
  ComplexD gamma, eta;  

  // r = r0 = chi - A psi_0
  A(r0, psi, isign);
  r[sub] = chi - r0;
  r0[sub] = r;

  flopcount.addFlops(A.nFlops());
  flopcount.addSiteFlops(2*Nc*Ns,sub);

 
  // f0 = A^\dag r_0 and plays a role like r_0 in inner products
  if(isign == PLUS) {
    A(f0,r0,MINUS);
  }
  else {
    A(f0,r0, PLUS);
  } 
  flopcount.addFlops(A.nFlops());

  // Now initialise v_0 = q_0 = z_0 = 0
  v[sub] = zero;
  q[sub] = zero;   
  z[sub] = zero;

  // u_0 = A r_0
  A(u,r,isign);   
  flopcount.addFlops(A.nFlops());


  // rho_0 := alpha_0 := omega_0 = 1
  rhon_1 = Double(1);
  alphan_1 = Double(1);
  omega = Double(1);

  // tau_0 = <r0,v0> = 0 since v0 = 0
  taun_1 = Double(0);

  // sigma_{-1} = 0 
  sigman_2 = Double(0);

  // pi_0 = <r0,q_0> = 0
  pi_c = Double(0);

  // sigma_{0} = <r0, Ar_0> = <r0, u_0>
  sigman_1 = innerProduct(r0,u,sub); 

  //**** NB: Paper says phi_0 = 0. 
  //**** PETSc code claims it is <r0,r0>
  //**** Using phi_0 = 0 leads to rho_1 = 0 
  //**** which leads to breakdown in iteration 2
  //**** So I will use the PETSc version which appears to work

  // phi_0 = <r0,r0>
  phi = innerProduct(r0,r0,sub);
  flopcount.addSiteFlops(16*Nc*Ns,sub);      // (for sigma & phi)


  // OK The big for loop
  for(int n = 1; n <= MaxBiCGStab && !convP ; n++) { 
 
    // Regular BiCGStab: rho = <r0,r> 
    // For IBiCGStab this has been unwound into 
    // the recurrencce:
    //
    //   rho = phi_{n-1} - omega_{n-1}*( sigma_{n-2} - alpha_{n-1}*pi_{n-1} )
    //
    // Check for rho=0. If it is it will lead to breakdown in next iteration
    rhon = phi - omega*(sigman_2 - alphan_1*pi_c);         // 16 flops


    if( toBool( real(rhon) == 0 ) && toBool( imag(rhon) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: rho = 0" << endl;
      QDP_abort(1);
    }
     
  
    // Regular BiCGStab: beta = ( rho_{n}/rho_{n-1})(alpha_{n-1}/omega_{n-1})
    // 
    // For IBiCGStab where one can use delta_n = beta*omega_{n-1}
    // it is useful to compute:
    //
    // delta_n = (rho_{n}/rho_{n-1})*alpha_{n-1}
    // beta_n  = delta_n/ omega_{n-1}
    ComplexD beta;
    ComplexD delta;
    delta =( rhon / rhon_1 ) * alphan_1;             // 15 flops
    beta = delta/omega;                              //  9 flops 


    // tau_n = <r0, v> needed for denominator of alpha 
    // but can be updated by recurrance
    taun = sigman_1 + beta*(taun_1- omega*pi_c);     // 16 flops 

    if( toBool( real(taun) == 0 ) && toBool( imag(taun) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: <r_0|v> = 0" << endl;
      QDP_abort(1);
    }

    // form alpha = rho/tau
    alphan = rhon / taun;                            // 9 flops


    // z_n plays role of alpha_n p_n in normal BiCGstab
    // it is only used to update the solution.
    //
    // NB one needs alpha_n p_{-1}
    //             = (alpha_n/alpha_{n-1}) (alpha_{n-1} p_{n-1})
    //             = (alpha_n/alpha_{n-1}) z_{n-1}
    // 
    // The Paper in line (12) of the algorithm leaves out this 
    // (alpha_n/alpha_{n-1}) factor.
    //
    // Also z update needs to be pulled before the v update (line 8) of paper
    // otherwise a shadow copy of v_{n-1} needs to be kept,
    //
    // z = alphan r_n-1 + (beta*alphan/alpha_{n-1})*z_{n-1}
    //                  - (beta*alphan*omegan_1)*v_{n-1}
    ComplexD bar = beta*alphan/alphan_1;           // 15 flops
    ComplexD alphadelta = alphan*delta;            //  6 flops 
    
#if 0
    tmp[sub] = bar*z;                           
    z[sub] = alphan*r+tmp;
    z[sub] -= alphadelta*v;                        // 22 Nc*Ns flops/site

    
    // v = u_{n-1} + beta*v_{n-1} - beta*omegan_1*q_n_1
    tmp[sub] = beta*v;  // 6Nc Ns
    v[sub] = u + tmp;   // 2Nc Ns
    v[sub] -= delta*q;  // 8Nc Ns
#else
    ibicgstab_zvupdates(r,z,v,u,q,alphan, bar, alphadelta, beta, delta, sub);
#endif


    // q = Av
    A(q,v,isign);

#if 0
    // t = u - alpha q
    t[sub] = u - alphan * q;         // 8 Nc Ns

    // s = r - alpha v
    s[sub] = r -  alphan*v;          // 8 Nc Ns


    // This should all be done with one sync point
    // BIG ALLREDUCE
    phi = innerProduct(r0,s,sub);         // 8 Nc Ns flops/site
    gamma = innerProduct(f0,s,sub);       // 8 Nc Ns flops/site
    pi_c = innerProduct(r0,q,sub);        // 8 Nc Ns flops/site
    eta = innerProduct(f0,t,sub);         // 8 Nc Ns flops/site
    theta = innerProduct(t,s,sub);        // 8 Nc Ns flops/site
    kappa = norm2(t,sub);                 // 4 Nc Ns flops/site
    rnorm_prev = norm2(r,sub);            // 4 Nc Ns flops/site 
#else
    ibicgstab_stupdates_reduces( alphan, 
				 r,
				 u,
				 v,
				 q,
				 r0,
				 f0,
				 s,
				 t,
				 phi,
				 pi_c, 
				 gamma,
				 eta,
				 theta, 
				 kappa,
				 rnorm_prev,
				 sub);
#endif

    // Collected flopcounts
    // coefficient recurrences:   flopcount.addFlops(86);
    // z & v updates              flopcount.addSiteFlops(38*Nc*Ns, sub);
    // q = Av                     flopcount.addFlops(A.nFlops());
    // s & t updates:             flopcount.addSiteFlops(16*Nc*Ns, sub);
    // 5 inner products, 2 norms: flopcount.addSiteFlops(48*Nc*Ns, sub);
    flopcount.addFlops(A.nFlops() + 86);
    flopcount.addSiteFlops(102*Nc*Ns, sub);

    // Check Norm: NB this is a lagged norm.
    // So we're checking whether we converged last iteration
    // This is the price of having only 1 allreduce
    if( toBool(rnorm_prev < rsd_sq ) ) {
      // Yes we've converged
      convP = true;
      ret.resid = sqrt(rnorm_prev);
      ret.n_count = n;

    }
    else { 

      // No we haven't converged ... continue...
      convP = false;

      // Check kappa for breakdown 
      if( toBool(kappa == 0) ) { 
	QDPIO::cerr << "Breakdown || Ms || = || t || = 0 " << endl;
	QDP_abort(1);
      }
      
      // Regular BiCGStab omega_n = <t,s> / <t,t> = theta/kappa
      omega = theta/kappa;                                     // 9 flops

      // sigma_n = <r0, A u_n> = gamma_n - omega_n * eta_n;
      //
      // NB: sigma_n is never explicitly used only sigma_{n-1}, sigma_{n-2}
      // So if I stuck this in sigma_n, I'd just end up moving it to sigma_{n-1}
      // So I'll just stick it straight into sigma_{n-1}
      sigman_2 = sigman_1;   // Preserve
      sigman_1 = gamma - omega*eta;                           // 8 flops


#if 0
      // Update r, psi
      // r = s - omega t
      r[sub] = s - omega*t;

      // psi = psi + omega s + z
      tmp[sub] = psi + omega*s;   
      psi[sub] = tmp + z;
#else
      ibicgstab_rxupdate(omega,s,t,z,r,psi,sub);
#endif

      // Recompute next u = A r
      A(u,r,isign); 


      // Update past values: Some of this could be saved I am sure
      rhon_1 = rhon;
      alphan_1 = alphan;
      taun_1 = taun;

      // Collected Flops
      // Omega + Sigma Updates: flopcount.addFlops(17);
      // r & x updates:         flopcount.addSiteFlops(18*Nc*Ns, sub)
      // u update               flopcount.addFlops(A.NFlops)
      flopcount.addFlops(A.nFlops()+17);
      flopcount.addSiteFlops(18*Nc*Ns,sub);
    }


  }
  
  swatch.stop();

  QDPIO::cout << "InvIBiCGStab: n = " << ret.n_count << " resid = " << ret.resid << endl;
  flopcount.report("invibicgstab", swatch.getTimeInSeconds());

  if ( ret.n_count == MaxBiCGStab ) { 
    QDPIO::cerr << "Nonconvergence of IBiCGStab. MaxIters reached " << endl;
  }

  finishKernels();
  return ret;
}


template<>
SystemSolverResults_t
InvIBiCGStab(const LinearOperator<LatticeFermionF>& A,
	    const LatticeFermionF& chi,
	    LatticeFermionF& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvIBiCGStab_a<LatticeFermionF, ComplexF>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

template<>
SystemSolverResults_t
InvIBiCGStab(const LinearOperator<LatticeFermionD>& A,
	    const LatticeFermionD& chi,
	    LatticeFermionD& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvIBiCGStab_a<LatticeFermionD, ComplexD>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

// Staggered
template<>
SystemSolverResults_t
InvIBiCGStab(const LinearOperator<LatticeStaggeredFermion>& A,
	    const LatticeStaggeredFermion& chi,
	    LatticeStaggeredFermion& psi,
	    const Real& RsdBiCGStab, 
	    int MaxBiCGStab, 
	    enum PlusMinus isign)

{
  return InvIBiCGStab_a<LatticeStaggeredFermion, Complex>(A, chi, psi, RsdBiCGStab, MaxBiCGStab, isign);
}

}  // end namespace Chroma
