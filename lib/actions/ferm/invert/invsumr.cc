#include "chromabase.h"
#include "actions/ferm/invert/invsumr.h"


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
void InvSUMR(const LinearOperator<T>& U,
	     const T& b,
	     T& x,
	     const Complex& zeta,
	     const Real& rho,
	     const Real& epsilon, 
	     int MaxSUMR, 
	     int& n_count)
{

  // First things first. We need r_0 = b - A x
  //
  // b - A x = b - (zeta I + rho U)x = b - zeta x - rho U x 
  //
  LatticeFermion r = b - zeta*x;
  LatticeFermion u;
  U(u, x, PLUS);
  r -= rho*u;


  // delta_0 = || r || 
  Real delta = sqrt(norm2(r));

  // phi_hat_1 = 1 / delta_0 
  Complex phi_hat = Real(1)/delta;

  // tau_hat = delta_0 / rho
  Complex tau_hat = delta/rho;

  LatticeFermion w, w_minus, p, nu_minus;

  // w_-1 = p_-1 = nu_0 = 0
  //
  // Iteration starts at m=1, so w_-1, p_-1 are in fact
  // w_{m-2}, p_{m-2}, nu_0 = nu_{m-1}
  w_minus=zero;
  p=zero;
  nu_minus = zero;

  
  Complex phi, s, lambda, r_m_minus;
  phi = Real(0);
  s   = Real(0);
  lambda = Real(0);
  r_m_minus = Real(0);

  Complex r_m_m, gamma;
  Real sigma;
  Complex c;

  r_m_m_prev = Real(1);
  gamma = Real(1);
  sigma = Real(1);
  c     = Real(1);

  LatticeFermion nu, nu_tilde;
  
  nu = r;
  nu_tilde = r;

  nu /= delta;
  nu_tilde /= delta;

  bool convP = false;

  Complex alpha;

  for(int m = 1; m <= MaxSUMR && !convP; m++) {

    // u := U nu_m
    U(u, nu, PLUS);

    // gamma_m := - nu_tilde*_m u = - nu_dagger u = - < nu, u>
    gamma = - innerProduct(nu_tilde, u);

    // sigma_m := ( ( 1 - | gamma_m | )(1 + | gamma_m |) )^{1/2}
    Real absgamma =sqrt( real(conj(gamma)*gamma ));
    sigma = sqrt( (Real(1) - absgamma)*(Real(1) + absgamma) );

    // alpha_m := - gamma_m  delta_m-1
    alpha = - gamma * delta;


    // r_m-1, m:= alpha_m * phi_hat_m + s_m-1*zeta/rho;
 
    r_m_minus = alpha * phi_hat + s*zeta/rho;

    Complex c_bar = conj(c);

    Complex r_hat_mm;

    r_hat_mm = alpha * phi_hat + c_bar*zeta/rho;

    Real abs_rhat_mm, abs_sigma_m, tmp_length;

    // | r_hat_mm |^2
    abs_rhat_mm = real(conj(r_hat_mm)*r_hat_mm);
    abs_sigma_m = sigma*sigma;
    tmp_length = sqrt(abs_rhat_mm + abs_sigma_m);

    c_hat = r_hat_mm / tmp_length;
    c = conj(c_hat);

    s = - sigma / tmp_length;

    r_m_m = -c*r_hat_mm + s*sigma;

    Complex tau = -c*tau_hat;
    tau_hat = s * tau_hat;

    Complex eta = tau / r_m_m;
    Complex kappa = r_m_minus/r_m_prev;
    
    r_m_prev = r_m_m;

    
    LatticeFermion w_minus_nu;
    w_minus_nu = w_minus - nu_minus;
    
    w_minus = w;
    w = alpha*p - kappa*w_minus_v;
    
    // p_{m-1} = p_{m-2} 
    p = p - w_minus_nu*lambda;

    w_minus_nu = w - nu;
    
    x = x - w_minus_nu*eta;

    if( sigma < fuzz ) {
      convP = true;
      continue;
    }

    delta = delta*sigma;

    phi = -c*phi_hat + s*conj(gamma)/delta;
    lambda = phi/r_m_m;
    phi_hat  = s* phi_hat + c_hat * conj(gamma)*delta;

    nu = (u + gamma*nu_tilde)/sigma;
    nu_tilde = (sigma*nu_tilde) + conj(gamma)*nu;

    Real taumod = real(conj(tau)*tau);
    if ( taumod < epsilon ) convP = true;

  }



}
