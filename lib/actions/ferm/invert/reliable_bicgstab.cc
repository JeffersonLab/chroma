// $Id: reliable_bicgstab.cc,v 3.9 2009-07-06 19:02:34 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/reliable_bicgstab.h"

#include "actions/ferm/invert/bicgstab_kernels.h"

namespace Chroma {

  using namespace BiCGStabKernels;

  template<typename T, typename TF, typename CF>
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

  BiCGStabKernels::initKernels();

  const Subset& s = A.subset();

  bool convP = false;

  // These are all the vectors. There should be
  // None declared later on. These declarations do 'mallocs'
  // under the hood. Want those out of the main loop.
  T b; 
  T tmp;
  T r_dble;
  T x_dble;

  TF r;
  TF r0;
  TF x; 
  TF p;
  TF v;
  TF t;

  int k;

  StopWatch swatch;
  FlopCounter flopcount;
  flopcount.reset();
  swatch.reset();
  swatch.start();

  x[s]=zero;
  p[s] = zero;
  v[s] = zero;

  Double rsd_sq =  Double(RsdBiCGStab)*Double(RsdBiCGStab)*norm2(chi,s);
  Double b_sq;


  A(tmp, psi, isign);

  // We could do all this in a onner
  // b_sq = minusTmpB(tmp, b, r, r0,s)
  //
  //b[s] = chi-tmp;
  //b_sq = norm2(b,s);
  // r[s] = b;

  xymz_normx(b,chi,tmp,b_sq,s);
  r[s] = b;
  r0[s] = b;
  Double r_sq = b_sq;
  QDPIO::cout << "r0 = " << b_sq << endl;;

  flopcount.addFlops(A.nFlops());
  flopcount.addSiteFlops(2*Nc*Ns,s);
  flopcount.addSiteFlops(4*Nc*Ns,s);  

  Double rNorm = sqrt(r_sq);
  Double r0Norm = rNorm;
  Double maxrx = rNorm;
  Double maxrr = rNorm;
  bool updateR = false;
  bool updateX = false;
  int rupdates = 0;
  int xupdates = 0;


  DComplex rho, rho_prev, alpha, omega;

  DComplex ctmp;
  Double t_norm;

  CF rho_r, alpha_r, omega_r;
  // rho_0 := alpha := omega = 1
  // Iterations start at k=1, so rho_0 is in rho_prev
  rho = Double(1);
  rho_prev = Double(1);
  alpha = Double(1);
  omega = Double(1);

  // The iterations 
  for(k = 0; k < MaxBiCGStab && !convP ; k++) { 
    
    if( k == 0 ) { 
      // I know that r_0 = r so <r_0|r>=norm2(r) = r_sq
      // rho = innerProduct(r0,r,s);
      rho = r_sq;
      p[s] = r;
    }
    else { 
      DComplex beta =(rho / rho_prev) * (alpha/omega);
      CF beta_r = beta;
      omega_r = omega;

      // NB: This could be done in a onner
      // rPlusBetaPMinusBetaOmegav(p, r, v, beta, omega, s) 

      // p = r + beta(p - omega v)
      // first work out p - omega v 
      // into tmp
      // then do p = r + beta tmp


      // tmp[s] = p - omega_r*v;
      // p[s] = r + beta_r*tmp;
      yxpaymabz(r, p, v, beta_r, omega_r, s);
      
    }

    // v = Ap
    AF(v,p,isign);

    // alpha = rho_{k+1} / < r_0 | v >
    // put <r_0 | v > into tmp
    ctmp = innerProduct(r0,v,s);

    if( toBool( real(ctmp) == 0 ) && toBool( imag(ctmp) == 0 ) ) {
      QDPIO::cout << "BiCGStab breakdown: <r_0|v> = 0" << endl;
      QDP_abort(1);
    }
    
    alpha = rho / ctmp;

    // Done with rho now, so save it into rho_prev
    rho_prev = rho;

    // s = r - alpha v
    // I can overlap s with r, because I recompute it at the end.
    alpha_r = alpha;
    // r[s]  -=  alpha_r*v;
    cxmay(r,v,alpha_r,s);


    // t = As  = Ar 
    AF(t,r,isign);


    // omega = < t | s > / < t | t > = < t | r > / norm2(t);
    // accumulate <t | s > = <t | r> into omega

    // As Mike tells me, I could do these together.
    // I can probably reduce these to a single ALLREDUCE/QMP_sum_double_array()
    // 
    // some routine like:  t_norm = normXCdotXY(t,r,s, iprod_r, iprod_i)
    // Double t_norm = norm2(t,s);
    // omega = innerProduct(t,r,s);

    
    norm2x_cdotxy(t,r, t_norm, omega, s);
    
    omega /= t_norm;

    // again
    // This is a simple xPlusAYPlusBz(x,r,p,omega,alpha)
    // psi = psi + omega s + alpha p 
    //     = psi + omega r + alpha p
    //
    // use tmp to compute psi + omega r
    // then add in the alpha p
    omega_r = omega;
    // tmp[s] = x + omega_r*r;
    // x[s] = tmp + alpha_r*p;


    xpaypbz(x,r,p,omega_r, alpha_r,s);

    // r = s - omega t = r - omega t1G

    // I can roll this all together 
    // r_sq = XMinusAYNormXCDotZX(r,t,r0,omega_r, omega_i, rho_r, rho_i, s),
    // r[s] -= omega_r*t;
    // r_sq = norm2(r,s);
    // rho = innerProduct(r0,r,s);

    xmay_normx_cdotzx(r, t, r0, omega_r, r_sq, rho,s); 

    // Flops so far: Standard BiCGStab Flops
    // -----------------------------------------
    flopcount.addSiteFlops(80*Nc*Ns,s);
    flopcount.addFlops(2*A.nFlops());
    // ------------------------------------------

    rNorm = sqrt(r_sq);

    if( toBool( rNorm > maxrx) ) maxrx = rNorm;
    if( toBool( rNorm > maxrr) ) maxrr = rNorm;

    updateX = toBool ( rNorm < Delta*r0Norm && r0Norm <= maxrx );
    updateR = toBool ( rNorm < Delta*maxrr && r0Norm <= maxrr ) || updateX;

    if( updateR ) { 
      // QDPIO::cout << "Iter " << k << ": updating r " << endl;
      rupdates++;
    
      x_dble[s] = x;

      A(tmp, x_dble, isign); // Use full solution so far

      // Roll this together - can eliminate r_dble which is an intermediary
	
      // r_dble[s] = b - tmp2;
      // r_sq = norm2(r_dble,s);
      // r[s] = r_dble;     
      xymz_normx(r_dble, b,tmp, r_sq,s);
      r[s]=r_dble;

      flopcount.addSiteFlops(6*Nc*Ns,s);
      flopcount.addFlops(A.nFlops());

      rNorm = sqrt(r_sq);
      maxrr = rNorm;
	
      
      if( updateX ) { 
	xupdates++;
	//QDPIO::cout << "Iter " << k << ": updating x " << endl;
	if( ! updateR ) { x_dble[s]=x; } // if updateR then this is done already
	psi[s] += x_dble; // Add on group accumulated solution in y
	flopcount.addSiteFlops(2*Nc*Ns,s);
	
	x[s] = zero; // zero y
	b[s] = r_dble;
	r0Norm = rNorm;
	maxrx = rNorm;
      }
      
    }


    if( toBool(r_sq < rsd_sq ) ) {
      
      convP = true;

      // if updateX true, then we have just updated psi
      // strictly x[s] should be zero, so it should be OK to add it
      // but why do the work if you don't need to
      x_dble[s] = x;
      psi[s]+=x_dble;
      flopcount.addSiteFlops(2*Nc*Ns,s);
      ret.resid = rNorm;
      ret.n_count = k;
    }
    else { 
      convP = false;
    }
    


  }
  swatch.stop();
  if( k >= MaxBiCGStab ) {
    QDPIO::cerr << "Nonconvergence of reliable BiCGStab. MaxIters = " << MaxBiCGStab << " exceeded" << endl;
    QDP_abort(1);
  }
  else { 
    QDPIO::cout << "reliable_bicgstab: n_count " << ret.n_count << " r-updates: " << rupdates << " xr-updates: " << xupdates  << endl;
    flopcount.report("reliable_bicgstab", swatch.getTimeInSeconds());
  }

  BiCGStabKernels::finishKernels();
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
  return RelInvBiCGStab_a<LatticeFermionF,LatticeFermionF, ComplexF>(A,A, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
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
  return RelInvBiCGStab_a<LatticeFermionD, LatticeFermionD, ComplexD>(A,A, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
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
  return RelInvBiCGStab_a<LatticeFermionD, LatticeFermionF, ComplexF>(A,AF, chi, psi, RsdBiCGStab, Delta, MaxBiCGStab, isign);
}


#if 0 

#endif 

}  // end namespace Chroma
