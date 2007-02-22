

#include "chromabase.h"
#include "actions/ferm/invert/inv_minres_array.h"

namespace Chroma {

template<typename T>
void InvMINRES_a(const LinearOperatorArray<T>& A,
		 const multi1d<T>& chi,
		 multi1d<T>& psi,
		 const Real& RsdCG,
		 int MaxCG,
		 int& n_count)
{
  START_CODE();

  // Compute v_1 = b - A x_0 
  int N5 = A.size();
  const Subset& s = A.subset();

  multi1d<T> v(N5);
  multi1d<T> Av(N5);

  A(Av, psi, PLUS);
  
  for(int n=0; n < N5; n++) { 
    v[n][s] = chi[n] - Av[n];
  }

  // beta_1 = || v_1 || ; eta = beta_1
  Real beta = sqrt(norm2(v,s));

  // r_norm_i = || b - A x || = beta_1
  Real r_norm = beta;
  Real b_norm = sqrt(norm2(chi, s));

  Real eta = beta;

  // gamma_1 = gamma_0 = 1
  Complex gamma, gamma_prev;
  gamma = Real(1); gamma_prev = Real(1);
  
  // sigma1 = sigma_0 = 0;
  Complex sigma , sigma_prev;
  sigma = Real(0); 
  sigma_prev = Real(0);

  // v_0 = 0; w_0 = w_1 = 0
  multi1d<T> v_prev(N5);

  multi1d<T> w_prev(N5);
  multi1d<T> w_2prev(N5);

  for(int n = 0; n < N5; n++) { 
    v_prev[n][s] = zero;
    w_prev[n][s] = zero;
    w_2prev[n][s] = zero;
  }

  bool convP = false;

  // Do the iteration proper
  for(int k = 1; k < MaxCG && !convP; k++) {

    // v_i = 1 / beta_i 
    Real ffactor = Real(1)/beta;
    for(int n=0; n < N5; n++) {
      v[n][s] *= ffactor;
    }

    // alpha_i = < v_i, A v_i >
    A(Av, v, PLUS);

    Complex alpha = Real(0);
    for(int n=0; n < N5; n++) { 
      alpha += innerProduct(v[n], Av[n], s);
    }

    
    // v_{i+1} = Av_{i} - alpha_i v_i - beta_i * v_{i-1}
    multi1d<T> v_plus(N5);
    for(int n=0; n < N5; n++) {
      v_plus[n][s] = Av[n] - alpha*v[n];
      v_plus[n][s] -= beta * v_prev[n];
    }
    
    // beta_{i+1} = || v_{i+1} ||
    // proper 5D norm I think
    Real beta_plus = sqrt(norm2(v_plus, s));

    
    // QR Part
    Complex delta;
    
    // delta = gamma_i alpha_i - gamma_{i-1}*sigma_{i}*beta_{i}
    delta = gamma * alpha - gamma_prev*sigma*beta;

    // rho_1 = sqrt( delta*2 + beta^2_{i+1}
    Real rho_1;
    rho_1 = sqrt( real(conj(delta)*delta) + beta_plus*beta_plus );

    // rho_2 = sigma_i*alpha_i + gamma_{i-1}*gamma_i*beta_i
    Complex rho_2;
    rho_2 = sigma*alpha + gamma_prev*gamma*beta;
    
    // rho_3 = sigma_{i-1} * beta_i
    Complex rho_3;
    rho_3 = sigma_prev * beta;

    Complex gamma_plus = delta / rho_1;
    Real sigma_plus = beta_plus / rho_1;

    multi1d<T> w(N5);
    // w_i = (v_i - rho_3*w_{i-2} - rho_2*w_{i-1})/rho_1
    for(int n=0; n < N5; n++) { 
      w[n][s]  = v[n] - rho_3*w_2prev[n];
      w[n][s] -= rho_2*w_prev[n];
      ffactor = Real(1)/rho_1;
      w[n][s] *= ffactor;
    }

    // x_{i} = x_{i-1} + gamma_{i+1}*eta* w_[i}
    // x in this case is psi
    Complex factor = gamma_plus*eta;
    for(int n=0; n < N5; n++) { 
      psi[n][s] += factor*w[n];
    }

    // || r_i || = | sigma_plus | || r_{i-1} ||
    r_norm *= fabs(sigma_plus) ;


    // Check convergence
    //  QDPIO::cout << "MINRES: iter " << k << " || r || = " << r_norm << endl;
    if( toBool( r_norm < RsdCG*b_norm ) ) {
      convP = true;
      END_CODE();
      return;
    }

    // eta = -sigma_{i+1}*eta
    eta *= -sigma_plus;

    // Now back everything up
    // v_prev = v
    // w_2prev = w_prev
    // w_prev = w
    //
    for(int n=0; n < N5; n++) { 
      v_prev[n][s] = v[n];
      v[n][s] = v_plus[n];
      w_2prev[n][s] = w_prev[n];
      w_prev[n][s] = w[n];
    }
 
    // beta_i <- beta_{i+1}
    beta = beta_plus;

    // gamma_{i-1} <- gamma_{i}
    gamma_prev = gamma;

    // gamma_{i} <- gamma_{i+1}
    gamma = gamma_plus;

    // sigma_{i-1} <- sigma_{i}
    sigma_prev = sigma;

    // sigma_{i} <- sigma_{i+1}
    sigma = sigma_plus;

    n_count = k;
  }

  END_CODE();
}


template<>
void InvMINRES(const LinearOperatorArray<LatticeFermion>& A,
	       const multi1d<LatticeFermion>& chi,
	       multi1d<LatticeFermion>& psi,
	       const Real& RsdCG,
	       int MaxCG,
	       int& n_count)
{
  InvMINRES_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

}  // end namespace Chroma
