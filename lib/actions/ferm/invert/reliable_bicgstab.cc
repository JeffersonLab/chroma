// $Id: reliable_bicgstab.cc,v 3.2 2009-05-20 18:22:34 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/reliable_bicgstab.h"

namespace Chroma {

  template<typename T, typename TF>
SystemSolverResults_t
RelInvBiCGStab_a(const LinearOperator<T>& A,
	      const LinearOperator<TF>& AF,
	      const T& chi,
	      T& psi,
	      const Real& RsdBiCGStab,
	      const Real& Delta,
	      int MaxBiCGStab, 
	      enum PlusMinus isign)
  {
  SystemSolverResults_t ret;

  const Subset& s = A.subset();

  bool convP = false;

  // First get r = r0 = chi - A psi
  TF r;
  T b; 
  T r_dble;
  T x_dble;

  b[s] = chi;
  {
    T tmp;
    A(tmp, psi, isign);
    b[s] -= tmp;
  }
  TF x; x[s]=zero;
 
  // now work out r= chi - Apsi = chi - r0
  r[s] = b; 
  TF r0; r0[s] = b;

  Double b_sq = norm2(b,s);
  Double r_sq = b_sq;
  Double rsd_sq =  Double(RsdBiCGStab)*Double(RsdBiCGStab)*b_sq;

  Double rNorm = sqrt(r_sq);
  Double r0Norm = rNorm;
  Double maxrx = rNorm;
  Double maxrr = rNorm;
  bool updateR = false;
  bool updateX = false;


  // Now initialise v = p = 0
  TF p;
  TF v;

  p[s] = zero;
  v[s] = zero;

  TF tmp;
  TF t;


  DComplex rho, rho_prev, alpha, omega;

  // rho_0 := alpha := omega = 1
  // Iterations start at k=1, so rho_0 is in rho_prev
  rho = Double(1);
  rho_prev = Double(1);
  alpha = Double(1);
  omega = Double(1);

  // The iterations 
  for(int k = 0; k <= MaxBiCGStab && !convP ; k++) { 
    
    // rho_{k+1} = < r_0 | r >
    if( k == 0 ) { 
      // I know that r_0 = r so <r_0|r>=norm2(r) = r_sq
      // rho = innerProduct(r0,r,s);
      rho = r_sq;
      p[s] = r;
    }
    else { 
      DComplex beta =(rho / rho_prev) * (alpha/omega);
    
      // p = r + beta(p - omega v)

      // first work out p - omega v 
      // into tmp
      // then do p = r + beta tmp
      tmp[s] = p - omega*v;
      p[s] = r + beta*tmp;
    }

    // v = Ap
#if 0
    { 
      TF vf; 
      TF pf;
      pf[s] = p;
      AF(vf,pf,isign);
      v[s] = vf;
    }
#else
    AF(v,p,isign);
#endif

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
    r[s]  -=  alpha*v;
    
    // t = As  = Ar 
#if 0
    { 
      TF tf;
      TF rf;
      rf[s] = r;

      AF(tf,rf,isign);
    
      t[s]=tf;
    }
#else
    AF(t,r,isign);
#endif

    // omega = < t | s > / < t | t > = < t | r > / norm2(t);

    // This does the full 5D norm
    // accumulate <t | s > = <t | r> into omega
    Double t_norm = norm2(t,s);
    omega = innerProduct(t,r,s);
    omega /= t_norm;

    // psi = psi + omega s + alpha p 
    //     = psi + omega r + alpha p
    //
    // use tmp to compute psi + omega r
    // then add in the alpha p
    tmp[s] = x + omega*r;
    x[s] = tmp + alpha*p;


    // r = s - omega t = r - omega t1G
    r[s] -= omega*t;


    r_sq = norm2(r,s);
    rho = innerProduct(r0,r,s);
    rNorm = sqrt(r_sq);

    if( toBool( rNorm > maxrx) ) maxrx = rNorm;
    if( toBool( rNorm > maxrr) ) maxrr = rNorm;

    updateX = toBool ( rNorm < Delta*r0Norm && r0Norm <= maxrx );
    updateR = toBool ( rNorm < Delta*maxrr && r0Norm <= maxrr ) || updateX;

    if( updateR ) { 
      // QDPIO::cout << "Iter " << k << ": updating r " << endl;
      
      {
	T tmp2;
	x_dble[s] = x;

	A(tmp2, x_dble, isign); // Use full solution so far
	r_dble[s] = b - tmp2;

	r[s] = r_dble;     // new R = b - Ax
	r_sq = norm2(r_dble,s);
      }
      rNorm = sqrt(r_sq);
      maxrr = rNorm;
	
      
      if( updateX ) { 
	//QDPIO::cout << "Iter " << k << ": updating x " << endl;
	if( ! updateR ) { x_dble[s]=x; } // if updateR then this is done already
	psi[s] += x_dble; // Add on group accumulated solution in y
	x[s] = zero; // zero y
	b[s] = r_dble;
	r0Norm = rNorm;
	maxrx = rNorm;
      }

    }
    // QDPIO::cout << "r_sq = " << r_sq << endl;

    //    QDPIO::cout << "Iteration " << k << " : r = " << r_norm << endl;
    if( toBool(r_sq < rsd_sq ) ) {
      
      convP = true;

      // if updateX true, then we have just updated psi
      // strictly x[s] should be zero, so it should be OK to add it
      // but why do the work if you don't need to
      x_dble[s] = x;
      psi[s]+=x_dble;
      

      ret.resid = rNorm;
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





SystemSolverResults_t
InvBiCGStabReliable(const LinearOperator<LatticeFermionF>& A,
		    const LatticeFermionF& chi,
		    LatticeFermionF& psi,
		    const Real& RsdBiCGStab, 
		    const Real& Delta,
		    int MaxBiCGStab, 
		    enum PlusMinus isign)

{
  return RelInvBiCGStab_a<LatticeFermionF,LatticeFermionF>(A,A, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}

  // Pure double
SystemSolverResults_t
InvBiCGStabReliable(const LinearOperator<LatticeFermionD>& A,
		    const LatticeFermionD& chi,
		    LatticeFermionD& psi,
		    const Real& RsdBiCGStab, 
		    const Real& Delta,
		    int MaxBiCGStab, 
		    enum PlusMinus isign)

{
  return RelInvBiCGStab_a<LatticeFermionD, LatticeFermionD>(A,A, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}

  // single double
SystemSolverResults_t
InvBiCGStabReliable(const LinearOperator<LatticeFermionD>& A,
		    const LinearOperator<LatticeFermionF>& AF,
		    const LatticeFermionD& chi,
		    LatticeFermionD& psi,
		    const Real& RsdBiCGStab, 
		    const Real& Delta,
		    int MaxBiCGStab, 
		    enum PlusMinus isign)

{
  return RelInvBiCGStab_a<LatticeFermionD, LatticeFermionF>(A,AF, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}


#if 0 

#endif 

}  // end namespace Chroma
