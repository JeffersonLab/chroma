#include "chromabase.h"
#include "actions/ferm/invert/invsumr.h"

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
void InvSUMR_a(const LinearOperator<T>& U,
	     const T& b,
	     T& x,
	     const Complex& zeta,
	     const Real& rho,
	     const Real& epsilon, 
	     int MaxSUMR, 
	     int& n_count)
{
  START_CODE();

  /********************************************************************/
  /* Initialisation                                                   */
  /********************************************************************/

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

  // If already converged then exit now.
  if( toBool( delta < epsilon*sqrt(norm2(b)) ) ) 
  {
    END_CODE();
    return;
  }

  // phi_hat_1 = 1 / delta_0 
  Complex phihat = Real(1)/delta;

  // tau_hat = delta_0 / rho
  Complex tauhat = delta/rho;

  LatticeFermion w_old, p, v_old;
  // w_-1 = p_-1 = v_0 = 0
  //
  // Iteration starts at m=1, so w_-1, p_-1 are in fact
  // w_{m-2}, p_{m-2}, v_{m-1}
  w_old=zero;
  p=w_old;
  v_old = w_old;

  // phi_0 = s_0 = lambda_0 := r_{-1,0} = 0
  Complex phi = Real(0);
  Real s   = Real(0);
  Complex lambda = Real(0);
  Complex r0 = Real(0);

  // r_{0,0} := gamma_0 := sigma_0 := c_0 := 1;
  Complex r1_old  = Real(1);
  Complex gamma = Real(1);
  Real sigma = Real(1);
  Complex c     = Real(1);

  // v_1 := \tilde{v}_1 := r_0 / delta_ 0
  LatticeFermion v, vtilde;
  Real ftmp = Real(1)/delta;
  v = ftmp*r;
  vtilde = ftmp*r;

  /***********************************************************/
  /* Start the iteration                                     */
  /***********************************************************/  
  bool convP = false;
  for(int iter = 1; iter <= MaxSUMR && !convP; iter++) {

    // u := U v_m
    U(u, v, PLUS);

    // gamma_m = - < vtilde, u >
    gamma = - innerProduct(vtilde, u);

    // sigma = ( ( 1 + | gamma | )(1 - | gamma | ) )^{1/2}
    //
    // NB I multiplied out the inner brackets
    // sigma = ( 1 - | gamma |^2 )^{1/2}
    sigma = sqrt( Real(1) - real( conj(gamma)*gamma ) );

    // alpha_m := - gamma_m  delta_{m-1}
    Complex alpha = - gamma * delta;

    // r_{m-1, m}:= alpha_m * phi_{m-1} + s_{m-1}*zeta/rho;
    r0 = alpha * phi + s*zeta/rho;

    // rhat_{m,m} := alpha_m phihat_m + conj(c_{m-1})*zeta/rho;
    Complex r1hat = alpha*phihat+conj(c)*zeta/rho;

    
    // conj(c_m) := rhat_{m,m}/( | rhat_{m,m} |^2 + | sigma_m |^2 )^{1/2}
    // s_m       := -sigma_m / ( | rhat_{m,m} |^2 + | sigma_m |^2 )^{1/2}
    //
    // Note I precompute the denominator into tmp_length.
    // and I conjugate the expression for conj(c) to get c directly
    // This trick is from the Wuppertal Group's MATLAB implementation
    Real abs_rhat_sq  = real(conj(r1hat)*r1hat);
    Real abs_sigma_sq = sigma*sigma;
    Real tmp_length   = sqrt(abs_rhat_sq + abs_sigma_sq);

    c = conj(r1hat)/tmp_length;
    s = - sigma / tmp_length;

    // r_{m,m} := -c_m*rhat_{m,m} + s_m*sigma_m;
    Complex r1 = -c*r1hat + cmplx(s*sigma,0);

    // tau_m   := -c_m*tauhat_m
    Complex tau = -c*tauhat;

    // tauhat_{m+1} := s_m tauhat_m
    tauhat = s * tauhat;

    // eta_m := tau_m / r_{m,m}
    Complex eta = tau / r1;

    // kappa_{m-1} := r_{m-1,m}/r_{m-1, m-1}
    Complex kappa = r0/r1_old;

    // keep r1
    r1_old=r1;

    // w_{m-1} := alpha_m * p_{m-2} - ( w_{m-2} - v_{m-1} ) kappa_{m-1}
    // p_{m-1} := p_{m-2} - ( w_{m-2} - v_{m-1} )*lambda_{m-1}
    //
    // NB: I precompute v_{m-2} - v_{m-1} into w_minus_v
    LatticeFermion w_minus_v;
    w_minus_v = w_old - v_old;

    LatticeFermion w = alpha*p - kappa*w_minus_v;
    p = p - lambda*w_minus_v;

    // x_{m-1} := x_{m-1} - ( w_{m-1} - v_m ) eta
    //
    // NB I compute w_{m-1} - v_m into w_minus_v

    w_minus_v = w-v;
    x = x - eta*w_minus_v;

    // I need to keep w as w_old for the next iteration
    w_old = w;

    // At this point the paper writes:
    // 
    // if (sigma_m = 0) then we have converged
    // so here I just go on and set the convergence flag in the else clause
    //
    if( toBool( sqrt(real(conj(sigma)*sigma)) > epsilon) ) {

      // delta_m := delta_{m-1} sigma_m
      delta = delta*sigma;

      // phi_m := -c_m phihat_m + s_m conj(gamma_m)/delta_m
      phi = -c*phihat + s*conj(gamma)/delta;

      // lambda_m := phi_m / r_{m,m}
      lambda = phi/r1;

      // phihat_{m+1} := s_m phihat_{m} + conj(c_m)*conj(gamma_m)/delta_m
      phihat  = s* phihat + conj(c) * conj(gamma)/delta;

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
      Real taumod = sqrt(real(conj(tauhat)*tauhat));
      if ( toBool(taumod < epsilon) ) convP = true;

      QDPIO::cout << "Iter " << iter << ":  | tauhat |=" << taumod << endl;
     
    }
    else { 

      // Else clause of the if sigma_m == 0. If sigma < epsilon then converged
      convP = true;
    }

    // Count the number of iterations
    n_count = iter;
  }

  // And we are done, either converged or not...
  if( n_count == MaxSUMR && ! convP ) { 
    QDPIO::cout << "Solver Nonconvergence Warning " << endl;
  }

  END_CODE();
}

template<>
void InvSUMR(const LinearOperator<LatticeFermion>& U,
	     const LatticeFermion& b,
	     LatticeFermion& x,
	     const Complex& zeta,
	     const Real& rho,
	     const Real& epsilon, 
	     int MaxSUMR, 
	     int& n_count)
{
  InvSUMR_a(U, b, x, zeta, rho, epsilon, MaxSUMR, n_count);
}

}  // end namespace Chroma
