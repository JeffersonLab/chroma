#include "chromabase.h"
#include "actions/ferm/invert/minv_rel_sumr.h"

namespace Chroma {

// Solve a shifted unitary system
//
// A x = b
//
// Where X is of the form: A = zeta I + rho U
//
// rho > 0 and zeta are complex, and U is unitary
//
// We solve with the method described in:
//
// "A Fast Minimal Residual Algorithm for Shifted Unitary Matrices"
// by Carl F. Jagels, and Lothar Reichel
// Numerical Linear Algebra with Applications, Vol 1(6), 555-570(1994)
//
// This paper is referenced by and applied to the Overlap Dirac Operator
// by G. Arnold, N. Cundy, J. van den Eshof, A Frommer, S. Krieg, T. Lippert,
// K. Schaefer "Numerical Methods for the QCD Overlap Operator: II. 
// Optimal Krylov Subspace Methods" -- hep-lat/0311025
// which is where the name SUMR was coined.
//

    
template<typename T>
void MInvRelSUMR_a(const LinearOperator<T>& U,
		   const T& b,
		   multi1d<T>& x,
		   const multi1d<Complex>& zeta,
		   const multi1d<Real>& rho,
		   const multi1d<Real>& epsilon, 
		   int MaxSUMR, 
		   int& n_count)
{
  START_CODE();

  // Sanity check
  int numroot = x.size();  
  if( zeta.size() != numroot ) { 
    QDPIO::cerr << "zeta size:"<< zeta.size()
		<<" is different from x.size:"<<numroot << endl;

    QDP_abort(1);
  }

  if( rho.size() != numroot) { 
    QDPIO::cerr << "rho size:"<< rho.size()
		<<" is different from x.size:"<<numroot << endl;

  }

  if( epsilon.size() != numroot) {
    QDPIO::cerr << "epsilon size:"<< epsilon.size()
		<<" is different from x.size:"<<numroot << endl;

  }

  // Auxiliary condition for stopping is when sigma < epsilon
  // To keep things consistent this should be the smallest epsilon passed in
  Real epsilon_sigma = epsilon[0];
  for(int shift = 0; shift < numroot; shift++) {
    if ( toBool( epsilon[shift] < epsilon_sigma ) ) {
      epsilon_sigma = epsilon[shift];
    }
  }
  
  /********************************************************************/
  /* Initialisation                                                   */
  /********************************************************************/

  // First things first. We need r_0 = b - A x
  // for all systems.
  // Restriction here, is that the initial vector x is zero
  // and so r_0 = b
  for(int shift = 0; shift < numroot; shift++) {
    x[shift] = zero;
  }

  LatticeFermion r;
  r = b;

  // delta_0_shift = || r_shift || 
  Real delta = sqrt(norm2(r));
  
  multi1d<Complex> phihat(numroot);
  multi1d<Complex> tauhat(numroot);

  // This vector is part of the Arnoldi process and should
  // be shift independent. 
  for(int shift = 0; shift < numroot; shift++) { 
    // phi_hat_1 = 1 / delta_0 
    phihat[shift] = Real(1)/delta;

    // tau_hat = delta_0 / rho
    tauhat[shift] = delta/rho[shift];
  }

  // v_0 = zero 
  LatticeFermion v_old = zero;

  // These are used in updating the individual solutions
  // and so are multi massed
  multi1d<LatticeFermion> w_old(numroot);
  multi1d<LatticeFermion> p(numroot);

  // w_-1 = p_-1 
  //
  // Iteration starts at m=1, so w_-1, p_-1 are in fact
  // w_{m-2}, p_{m-2}, v_{m-1}
  for(int shift = 0; shift < numroot; shift++) { 
    w_old[shift]=zero;
    p[shift]=zero;
  }


  // These are part of the Arnoldi process and ar shift independent
  // gamma_0:= 1, sigma := 1
  Complex gamma = Real(1);
  Real sigma = Real(1);

  // These are all shift dependent
  multi1d<Complex> phi(numroot);
  multi1d<Real>    s(numroot);
  multi1d<Complex> lambda(numroot);
  multi1d<Complex> r0(numroot);
  multi1d<Complex> r1_old(numroot);
  multi1d<Complex> c(numroot);

  // phi_0 = s_0 = lambda_0 := r_{-1,0} = 0  
  // r_{0,0} := c_0 := 1;
  for(int shift = 0; shift < numroot; shift++) { 
    phi[shift] = Real(0);
    s[shift]   = Real(0);
    lambda[shift] = Real(0);
    r0[shift] = Real(0);
    r1_old[shift] = Real(1);
    c[shift] = Real(1);
  }
 
  // It is here, that we couple the individual systems. 
  // Because r_0 is the same for all the systems
  // v_1 is the same for all the system

  // v_1 := \tilde{v}_1 := r_0 / delta_ 0
  LatticeFermion v;
  LatticeFermion vtilde;

  Real ftmp = Real(1)/delta;
  v = ftmp*r;
  vtilde = ftmp*r;

  /***********************************************************/
  /* Start the iteration                                     */
  /***********************************************************/  
  multi1d<bool> convP(numroot);
  convP = false;
  bool allConvP = false;


  for(int iter = 1; iter <= MaxSUMR && !allConvP; iter++) {
    LatticeFermion u;

    // Updating u, gamma and sigma are common to all systems.

    // the inner solver precision is supposed to be 
    // epsilon / || r ||
    //
    // we get the bound for || r || from tauhat.
    // and here we need tauhat from the system with the smallest
    // unconverged shift
    int unc=0;
    while( convP[unc] == true && unc < numroot) { unc++; }
    
    // If there are other shifts find the smallest unconverged one.
    if( unc < numroot ) { 
      Real mod_me = sqrt(real(conj(zeta[unc])*zeta[unc]));

      // look through the other shifts
      for(int shift=unc+1; shift < numroot; shift++) { 

	// Only compare if they are  unconverged 
	if( convP[shift] == false ) {
	  Real mod_trial = sqrt(real(conj(zeta[shift])*zeta[shift]));
	  
	  if ( toBool(mod_trial < mod_me) ) { 
	    // The trial shift has smaller norm than me and is unconverged
	    // so we take it as the smallest and save its mod
	    unc = shift;
	    mod_me = mod_trial;
	  }

	}
      }
    }
    else {
      QDPIO::cerr << "All systems appear to be converged. I shouldnt be here" << endl;
    }
   
    // unc now contains the shift index of the system with the smallest
    // shift that is still unconverged. Use the tauhat of this for the
    // relaxation.

    Real inner_tol = epsilon_sigma / sqrt(real(conj(tauhat[unc])*tauhat[unc]));
   
      

    // u := U v_m
    U(u, v, PLUS, inner_tol);

    // gamma_m = - < vtilde, u >
    gamma = - innerProduct(vtilde, u);

    // sigma = ( ( 1 + | gamma | )(1 - | gamma | ) )^{1/2}
    //
    // NB I multiplied out the inner brackets
    // sigma = ( 1 - | gamma |^2 )^{1/2}
    sigma = sqrt( Real(1) - real( conj(gamma)*gamma ) );

    // multi1d<Complex> alpha(numroot);
    Complex alpha;
    alpha = -gamma*delta;
  
    multi1d<Complex> r1hat(numroot);

    // the abs_rhat_sq abs_sigma_sq and tmp_length are truly temporary
    // further, sigma is shift independent so I compute it here.
    Real abs_rhat_sq;
    Real abs_sigma_sq = sigma*sigma;
    Real tmp_length;
    multi1d<Complex> r1(numroot);

    // Go through all the shifts and update w, p, and x 
    for(int shift = 0; shift < numroot; shift++) {

      // Only do unconverged systems
      if( !convP[shift] ) {

	// r_{m-1, m}:= alpha_m * phi_{m-1} + s_{m-1}*zeta/rho;
	r0[shift] = alpha * phi[shift] + s[shift]*zeta[shift]/rho[shift];

	// rhat_{m,m} := alpha_m phihat_m + conj(c_{m-1})*zeta/rho;
	r1hat = alpha * phihat[shift]+conj(c[shift])*zeta[shift]/rho[shift];

	// conj(c_m) := rhat_{m,m}/( | rhat_{m,m} |^2 + | sigma_m |^2 )^{1/2}
	// s_m       := -sigma_m / ( | rhat_{m,m} |^2 + | sigma_m |^2 )^{1/2}
	//
	// Note I precompute the denominator into tmp_length.
	// and I conjugate the expression for conj(c) to get c directly
	// This trick is from the Wuppertal Group's MATLAB implementation
	abs_rhat_sq  = real(conj(r1hat[shift])*r1hat[shift]);
	tmp_length   = sqrt(abs_rhat_sq + abs_sigma_sq);

	c[shift] = conj(r1hat[shift])/tmp_length;
	s[shift] = - sigma / tmp_length;

	// r_{m,m} := -c_m*rhat_{m,m} + s_m*sigma_m;
	r1[shift] = -c[shift]*r1hat[shift] + cmplx(s[shift]*sigma,0);

	// tau_m   := -c_m*tauhat_m
	Complex tau = -c[shift]*tauhat[shift];

	// tauhat_{m+1} := s_m tauhat_m
	tauhat[shift] = s[shift] * tauhat[shift];

	// eta_m := tau_m / r_{m,m}
	Complex eta = tau / r1[shift];

	// kappa_{m-1} := r_{m-1,m}/r_{m-1, m-1}
	Complex kappa = r0[shift]/r1_old[shift];

	// keep r1
	r1_old[shift]=r1[shift];

	// w_{m-1} := alpha_m * p_{m-2} - ( w_{m-2} - v_{m-1} ) kappa_{m-1}
	// p_{m-1} := p_{m-2} - ( w_{m-2} - v_{m-1} )*lambda_{m-1}
	//
	// NB: I precompute v_{m-2} - v_{m-1} into w_minus_v
	LatticeFermion w_minus_v;
	w_minus_v = w_old[shift] - v_old;

	LatticeFermion w = alpha*p[shift] - kappa*w_minus_v;
	p[shift] = p[shift] - lambda[shift]*w_minus_v;

	// x_{m-1} := x_{m-1} - ( w_{m-1} - v_m ) eta
	//
	// NB I compute w_{m-1} - v_m into w_minus_v
	
	w_minus_v = w-v;
	x[shift] = x[shift] - eta*w_minus_v;

	// I need to keep w as w_old for the next iteration
	w_old[shift] = w;

      }
    }

    // At this point the paper writes:
    // 
    // if (sigma_m = 0) then we have converged
    // 
    // In this case I regard all systems as converted
    //
    if( toBool( sqrt(real(conj(sigma)*sigma)) > epsilon_sigma ) ) {

      // Here we haven't converged
      
      // delta_m := delta_{m-1} sigma_m
      delta = delta*sigma;
      
      // Update phi, lambda, phihat
      for(int shift=0; shift < numroot; shift++) { 

	// But only for the update systems
	if( !convP[shift] ) { 

	  // phi_m := -c_m phihat_m + s_m conj(gamma_m)/delta_m
	  phi[shift] = -c[shift]*phihat[shift]
	    + s[shift]*conj(gamma)/delta;
	  
	  // lambda_m := phi_m / r_{m,m}
	  lambda[shift] = phi[shift]/r1[shift];
	  
	  // phihat_{m+1} := s_m phihat_{m} + conj(c_m)*conj(gamma_m)/delta_m
	  phihat[shift]  = s[shift]* phihat[shift] 
	    + conj(c[shift]) * conj(gamma)/delta;
	  
	}
      }

      // preserve v as v_old
      v_old = v;
      
      // v_{m+1} := (1/sigma)( u + gamma_m vtilde_m );
      v = u + gamma*vtilde;
      Real ftmp2 = Real(1)/sigma;
      v *= ftmp2;

      // vtilde_{m+1} := sigma_m vtilde_m + conj(gamma_m)*v_{m+1} 
      Complex gconj = conj(gamma);
      Complex csigma=cmplx(sigma,0);
      vtilde = csigma*vtilde + gconj*v;

      // Normalise vtilde -- I found this in the Wuppertal MATLAB code
      // It is not prescribed by the Reichel/Jagels paper
      Real ftmp3 = Real(1)/sqrt(norm2(vtilde));
      vtilde *= ftmp3;

      // Check whether we have converged or not:
      //
      // converged if  | tauhat | < epsilon 
      //
      // Assume we have all converged, and we BOOLEAN AND
      // with individual systems
      allConvP = true;

      // Go through all the shifted systems
      for(int shift=0; shift < numroot; shift++) { 

	// Only check unconverged systems
	if( ! convP[shift] ) { 
	 
	  // Get | tauhat |
	  Real taumod = sqrt(real(conj(tauhat[shift])*tauhat[shift]));
	  
	  QDPIO::cout << "Iter " << iter << ":  Shift: " << shift 
		      <<" | tauhat |=" << taumod << endl;
	  
	  // Check convergence
	  if ( toBool(taumod < epsilon[shift] ) ) convP[shift] = true;

	  // Boolean AND with Overall convergence flag
	  allConvP &= convP[shift];
	}
	
      }
    }
    else { 
      
      // Else clause of the if sigma_m == 0. If sigma < epsilon then converged
      allConvP = true;
      
    }
    
    // Count the number of iterations
    n_count = iter;
  }
  
  // And we are done, either converged or not...
  if( n_count == MaxSUMR && ! allConvP ) { 
    QDPIO::cout << "Solver Nonconvergence Warning " << endl;
    QDP_abort(1);
  }

  END_CODE();
}

template<>
void MInvRelSUMR(const LinearOperator<LatticeFermion>& U,
		 const LatticeFermion& b,
		 multi1d<LatticeFermion>& x,
		 const multi1d<Complex>& zeta,
		 const multi1d<Real>& rho,
		 const multi1d<Real>& epsilon, 
		 int MaxSUMR, 
		 int& n_count)
{
  MInvRelSUMR_a(U, b, x, zeta, rho, epsilon, MaxSUMR, n_count);
}

}  // end namespace Chroma
