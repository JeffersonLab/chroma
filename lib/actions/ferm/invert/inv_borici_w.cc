#include "chromabase.h"
#include "actions/ferm/invert/inv_borici_w.h"
#include "actions/ferm/invert/inv_minres_array.h"

namespace Chroma {


template<typename T>
void InvBorici_a( const LinearOperator<T>& D_4,
		  const LinearOperatorArray<T>& D_5,
		  const LinearOperatorArray<T>& D_dag_D_5,
		  const T& b,
		  T& x,
		  const Real& tol,
		  const Real& tol_1,
		  const int MaxIter,
		  const int MaxIter5D,
		  const Real& m,
		  int& n_iters)
{
  START_CODE();

  // Initial setup
  // x_0 = 0, r_0 = b, tol_0 = 1
  x = zero;
  T r = b;
  
  // Other useful things
  Double b_norm = sqrt(norm2(b));
  Double r_norm;

  // Convergence flag
  bool convP = false;


  // The 5D Auxiliaries
  int N5 = D_5.size();
  multi1d<T> chi(N5);
  multi1d<T> eta(N5);
  //multi1d<T> tmp(N5);
  int n_iters5d;

  int G5 = Ns*Ns-1;
  T y;

  Real tol_0 = 1;

  for(int k = 1; k <= MaxIter && !convP; k++) { 


    tol_0 *= tol_1;

    // Now solve the 5D system to tol_1
    for(int n=0; n < N5-1; n++) {
      eta[n] = zero;
      //      tmp[n] = zero;
    }
    // RHS is r
    eta[N5-1] = Gamma(G5)*r;
    // tmp[N5-1] = zero;


   
    // Solve D_5 chi = eta 
    // or    D_5 y = r 
    // 
    // Solve with CGNE   D^{dag} D tmp = r
    //  then            chi = D^{dag} tmp = D^{dag} (D^{dag}D)^{-1} r
    //                                    = D^{-1} r

    InvMINRES(D_5, eta, chi, tol_0, MaxIter5D, n_iters5d);
    // D_5(chi,tmp,MINUS);
    
    y = Real(2)/(Real(1)-m)*chi[N5-1];

    // Keep y as next initial guess
    // x_i = x_{i-1} + y 
    x += y;
    
    // r_i = b - D_ov x_i
    T tmp4;
    D_4(tmp4, x, PLUS);

    r =b - tmp4;

    // Monitoring
    r_norm = sqrt(norm2(r));
    QDPIO::cout << "InvBorici: iter " << k 
		<< " || r || = " << r_norm 
		<< " || b || = " << b_norm 
		<< " || r ||/|| b || = " << r_norm / b_norm 
		<< " target = " << tol << endl;
    
    if( toBool( r_norm < tol*b_norm ) ) {
      convP = true;
    }
    n_iters = k;

  }

  if( !convP ) {
    QDPIO::cerr << "Nonconvergence warning: " << n_iters << " iters. || r || / || b || = " << r_norm / b_norm;
  }

  END_CODE();
}

template<>
void InvBorici( const LinearOperator<LatticeFermion>& D_4,
		const LinearOperatorArray<LatticeFermion>& D_5,
		const LinearOperatorArray<LatticeFermion>& D_dag_D_5,
		const LatticeFermion& b,
		LatticeFermion& x,
		const Real& tol,
		const Real& tol_1,
		const int MaxIter,
		const int MaxIter5D,
		const Real& m,
		int& n_iters)
{
  InvBorici_a(D_4, D_5, D_dag_D_5, b, x, tol, tol_1, MaxIter, MaxIter5D, m, n_iters);
}

}  // end namespace Chroma
