// $Id: reliable_ibicgstab.cc,v 3.2 2009-07-07 19:13:20 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/reliable_ibicgstab.h"

#include "actions/ferm/invert/bicgstab_kernels.h"

namespace Chroma {

  using namespace BiCGStabKernels;

  template<typename T, typename TF, typename CF>
SystemSolverResults_t
RelInvIBiCGStab_a(const LinearOperator<T>& A,
	      const LinearOperator<TF>& AF,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      const Real& Delta,
	      int MaxBiCGStab, 
	      enum PlusMinus isign)
  {
  SystemSolverResults_t ret;

  BiCGStabKernels::initKernels();
  StopWatch swatch;
  FlopCounter flopcount;
  flopcount.reset();
  const Subset& sub = A.subset();
  bool convP = false;

  T b; 
  T dtmp;
  T r_dble;
  T x_dble;


  TF r;
  TF r0;
  TF v;
  TF tmp;
  TF t;
  TF s;
  TF u;
  TF q;
  TF f0; 
  TF z;
  TF x;

  TF vn_1, zn_1,qn_1;

  int k;
  ComplexD rhon, rhon_1;             // rho_n, rho_{n-1}
  ComplexD alphan, alphan_1;         // alpha_{n}, alpha_{n-1}
  ComplexD omegan, omegan_1;          // omega_{n} AND omega_{n-1} (before omega update)
  ComplexD taun,taun_1;              // tau_{n}, tau_{n-1}
  RealD kappa;                       
  RealD rnorm_prev;                  // lagged residuum
  ComplexD theta;
  ComplexD sigman_2, sigman_1;       
  ComplexD phin, phin_1;                      // phi_{n} AND phi_{n-1} before phi update
  ComplexD pin, pin_1;                     // pi_{n} and pi_{n-1} before pi update
  ComplexD gamma, eta;  
  swatch.reset();
  swatch.start();


  Double rsd_sq =  RsdBiCGStab*RsdBiCGStab*norm2(chi,sub);
  Double b_sq;
  flopcount.addSiteFlops(4*Nc*Ns,sub);

  // Now initialise x_0 = v_0 = q_0 = z_0 = 0 
  x[sub] = zero;
  vn_1[sub] = zero;
  qn_1[sub] = zero;   
  zn_1[sub] = zero;



  // r = r0 = chi - A psi_0
  A(dtmp, psi, isign);
  xymz_normx(b,chi,dtmp,b_sq,sub);
  r[sub] = b;
  r0[sub] = b;

  QDPIO::cout << "r0 = " << b_sq << endl;;

  flopcount.addFlops(A.nFlops());
  flopcount.addSiteFlops(2*Nc*Ns,sub);
  flopcount.addSiteFlops(4*Nc*Ns,sub);  

  Double rNorm = sqrt(b_sq);
  Double r0Norm = rNorm;
  Double maxrx = rNorm;
  Double maxrr = rNorm;
  bool updateR = false;
  bool updateX = false;
  int xupdates = 0;
  int rupdates = 0;

  // f0 = A^\dag r_0 and plays a role like r_0 in inner products
  if(isign == PLUS) {
    AF(f0,r0,MINUS);
  }
  else {
    AF(f0,r0, PLUS);
  } 
  flopcount.addFlops(A.nFlops());


  // u_0 = A r_0
  AF(u,r,isign);   
  flopcount.addFlops(A.nFlops());


  // rho_0 := alpha_0 := omega_0 = 1
  rhon_1 = Double(1);
  alphan_1 = Double(1);
  omegan_1 = Double(1);

  // tau_0 = <r0,v0> = 0 since v0 = 0
  taun_1 = Double(0);

  // sigma_{-1} = 0 
  sigman_2 = Double(0);

  // pi_0 = <r0,q_0> = 0
  pin_1 = Double(0);

  // sigma_{0} = <r0, Ar_0> = <r0, u_0>
  sigman_1 = innerProduct(r0,u,sub); 

  //**** NB: Paper says phi_0 = 0. 
  //**** PETSc code claims it is <r0,r0>
  //**** Using phi_0 = 0 leads to rho_1 = 0 
  //**** which leads to breakdown in iteration 2
  //**** So I will use the PETSc version which appears to work

  // phi_0 = <r0,r0>
  phin_1 = innerProduct(r0,r0,sub);
  flopcount.addSiteFlops(16*Nc*Ns,sub);      // (for sigma & phi)

  /********** ITERATION STARTS HERE ***************/

  // OK The big for loop
  for(int n = 1; n <= MaxBiCGStab && !convP ; n++) { 
 
    // Regular BiCGStab: rho = <r0,r> 
    // For IBiCGStab this has been unwound into 
    // the recurrencce:
    //
    //   rho = phi_{n-1} - omega_{n-1}*( sigma_{n-2} - alpha_{n-1}*pi_{n-1} )
    //
    // Check for rho=0. If it is it will lead to breakdown in next iteration
    rhon = phin_1 - omegan_1*(sigman_2 - alphan_1*pin_1);         // 16 flops


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
    beta = delta/omegan_1;                              //  9 flops 


    // tau_n = <r0, v> needed for denominator of alpha 
    // but can be updated by recurrance
    taun = sigman_1 + beta*(taun_1- omegan_1*pin_1);     // 16 flops 

    if( toBool( real(taun) == 0 ) && toBool( imag(taun) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: n="<<n<<" <r_0|v> = 0" << endl;
      QDPIO::cout << "sigman_1 = " << sigman_1 << endl;
      QDPIO::cout << "beta= " << beta << endl;
      QDPIO::cout << "taun_1 = " << taun_1 << endl;
      QDPIO::cout << "ometan_1 = " << omegan_1 << endl;
      QDPIO::cout << "pin_1 = " << pin_1 << endl;

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
    tmp[sub] = bar*zn_1;                           
    z[sub] = alphan*r+tmp;
    z[sub] -= alphadelta*vn_1;                        // 22 Nc*Ns flops/site

    
    // v = u_{n-1} + beta*v_{n-1} - beta*omegan_1*q_n_1
    tmp[sub] = beta*vn_1;  // 6Nc Ns
    v[sub] = u + tmp;   // 2Nc Ns
    v[sub] -= delta*qn_1;  // 8Nc Ns
#else
    v[sub]=vn_1;
    z[sub]=zn_1;

    ibicgstab_zvupdates(r,z,v,u,qn_1,alphan, bar, alphadelta, beta, delta, sub);
#endif

    // q = Av
    AF(q,v,isign);

#if 0
    // t = u - alpha q
    t[sub] = u - alphan * q;         // 8 Nc Ns

    // s = r - alpha v
    s[sub] = r -  alphan*v;          // 8 Nc Ns


    // This should all be done with one sync point
    // BIG ALLREDUCE


    phin = innerProduct(r0,s,sub);         // 8 Nc Ns flops/site
    gamma = innerProduct(f0,s,sub);       // 8 Nc Ns flops/site
    pin = innerProduct(r0,q,sub);        // 8 Nc Ns flops/site
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
				 phin,
				 pin, 
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

 
    // Check Norm: See if we converged in the last iteration
    // If so go to exit. This is yucky but makes the rest 
    // of the logic easier
    if( toBool(rnorm_prev < rsd_sq ) ) {
      // Yes we've converged
      convP = true;

      // if updateX true, then we have just updated psi
      // strictly x[sub] should be zero, so it should be OK to add it
      // but why do the work if you don't need to
      x_dble[sub] = x;
      psi[sub]+=x_dble;
      flopcount.addSiteFlops(2*Nc*Ns,sub);
      ret.resid = sqrt(rnorm_prev);
      ret.n_count = n;
      goto exit;
    }
 

#if 1
    // Begin Reliable Updating ideas...
    rNorm = sqrt(rnorm_prev);
    if( toBool( rNorm > maxrx) ) maxrx = rNorm;
    if( toBool( rNorm > maxrr) ) maxrr = rNorm;
    
    updateX = toBool ( rNorm < Delta*r0Norm && r0Norm <= maxrx );
    updateR = toBool ( rNorm < Delta*maxrr && r0Norm <= maxrr ) || updateX;

    if( updateR ) { 

      // Replace last residuum
      x_dble[sub] = x;

      A(dtmp, x_dble, isign); // Use full solution so far

      // r_dble[sub] = b - tmp2;
      // r_sq = norm2(r_dble,sub);
      // r[s] = r_dble;     
      xymz_normx(r_dble, b,dtmp, rnorm_prev,sub);
      r[sub]=r_dble;
      flopcount.addFlops(A.nFlops());
      flopcount.addSiteFlops(6*Nc*Ns,sub);


      // Must also reset un_1
      AF(u,r,isign);
      
      // Recomputing sigma_{n-1} = < r0, u_{n-1} > 
      // from the new u_{n-1} really helps stability
      sigman_1 = innerProduct(r0,u,sub); 
      flopcount.addFlops(A.nFlops());
      flopcount.addSiteFlops(8*Nc*Ns,sub);


      rNorm = sqrt(rnorm_prev);
      maxrr = rNorm;
      rupdates++;
      
      if( updateX ) { 
	//QDPIO::cout << "Iter " << k << ": updating x " << endl;
	if( ! updateR ) { x_dble[sub]=x; } // if updateR then this is done already
	psi[sub] += x_dble; // Add on group accumulated solution in y
	flopcount.addSiteFlops(2*Nc*Ns,sub);
	
	x[sub] = zero; // zero y
	b[sub] = r_dble;
	r0Norm = rNorm;
	maxrx = rNorm;
	xupdates++;
      }
      
    }

#endif

    if( ! updateR ) {

      // Carry on with this iteration
      // Check kappa for breakdown 
      if( toBool(kappa == 0) ) { 
	QDPIO::cerr << "Breakdown || Ms || = || t || = 0 " << endl;
	QDP_abort(1);
      }
      
      // Regular BiCGStab omega_n = <t,s> / <t,t> = theta/kappa
      omegan = theta/kappa;                                     // 9 flops
      
      
      
#if 0
      // Update r, x
      // r = s - omega t
      r[sub] = s - omegan*t;
      
      // x = x + omega s + z
      tmp[sub] = x + omegan*s;   
      x[sub] = tmp + z;
#else
      ibicgstab_rxupdate(omegan,s,t,z,r,x,sub);
#endif
      
      // Recompute next u = A r
      AF(u,r,isign); 
      
      // sigma_n = <r0, A u_n> = gamma_n - omega_n * eta_n;
      //
      // NB: sigma_n is never explicitly used only sigma_{n-1}, sigma_{n-2}
      // So if I stuck this in sigma_n, I'd just end up moving it to sigma_{n-1}
      // So I'll just stick it straight into sigma_{n-1}
      sigman_2 = sigman_1;   // Preserve
      sigman_1 = gamma - omegan*eta;                           // 8 flops
      
      // Update past values: Some of this could be saved I am sure
      rhon_1 = rhon;
      alphan_1 = alphan;
      taun_1 = taun;
      omegan_1 = omegan;
      pin_1 = pin;
      phin_1 = phin;
     
      vn_1[sub]=v;
      qn_1[sub]=q;
      zn_1[sub]=z;
       
      // Collected Flops
      // Omega + Sigma Updates: flopcount.addFlops(17);
      // r & x updates:         flopcount.addSiteFlops(18*Nc*Ns, sub)
      // u update               flopcount.addFlops(A.NFlops)
      flopcount.addFlops(A.nFlops()+17);
      flopcount.addSiteFlops(18*Nc*Ns,sub);
    }
  }
  
exit:  swatch.stop();

  QDPIO::cout << "InvIBiCGStabReliable: n = " << ret.n_count << " r-updates: " << rupdates << " xr-updates: " << xupdates << endl;
 
  flopcount.report("reliable_invibicgstab", swatch.getTimeInSeconds());

  if ( ret.n_count == MaxBiCGStab ) { 
    QDPIO::cerr << "Nonconvergence of IBiCGStab. MaxIters reached " << endl;
  }


  BiCGStabKernels::finishKernels();
  return ret;

}





SystemSolverResults_t
InvIBiCGStabReliable(const LinearOperator<LatticeFermionF>& A,
		    const LatticeFermionF& chi,
		    LatticeFermionF& psi,
		    const Real& RsdBiCGStab, 
		    const Real& Delta,
		    int MaxBiCGStab, 
		    enum PlusMinus isign)

{
  return RelInvIBiCGStab_a<LatticeFermionF,LatticeFermionF, ComplexF>(A,A, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}

  // Pure double
SystemSolverResults_t
InvIBiCGStabReliable(const LinearOperator<LatticeFermionD>& A,
		    const LatticeFermionD& chi,
		    LatticeFermionD& psi,
		    const Real& RsdBiCGStab, 
		    const Real& Delta,
		    int MaxBiCGStab, 
		    enum PlusMinus isign)

{
  return RelInvIBiCGStab_a<LatticeFermionD, LatticeFermionD, ComplexD>(A,A, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}

  // single double
SystemSolverResults_t
InvIBiCGStabReliable(const LinearOperator<LatticeFermionD>& A,
		    const LinearOperator<LatticeFermionF>& AF,
		    const LatticeFermionD& chi,
		    LatticeFermionD& psi,
		    const Real& RsdBiCGStab, 
		    const Real& Delta,
		    int MaxBiCGStab, 
		    enum PlusMinus isign)

{
  return RelInvIBiCGStab_a<LatticeFermionD, LatticeFermionF, ComplexF>(A,AF, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}


#if 0 

#endif 

}  // end namespace Chroma
