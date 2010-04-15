// $Id: reliable_cg.cc,v 3.3 2009-06-01 16:24:54 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/reliable_cg.h"

namespace Chroma {

  template<typename T, typename TF, typename RF>
SystemSolverResults_t
RelInvCG_a(const LinearOperator<T>& A,
	   const LinearOperator<TF>& AF,
	   const T& chi,
	   T& psi,
	   const Real& RsdCG,
	   const Real& Delta,
	   int MaxCG)
  {
    START_CODE();
    SystemSolverResults_t ret;

    const Subset& s = A.subset();
    
    bool convP = false;

    // First get r = r0 = chi - A psi
    TF r;
    T b; 
    T r_dble;
    T x_dble;
    int k;

    StopWatch swatch;
    FlopCounter flopcount;
    flopcount.reset();
    swatch.reset();
    swatch.start();

    

    b[s] = chi;
    Double chi_norm = norm2(chi,s);
    Double rsd_sq=RsdCG*RsdCG*chi_norm;

    {
      T tmp1, tmp2;
      A(tmp1, psi, PLUS);
      A(tmp2, tmp1, MINUS);
      b[s] -= tmp2;
      flopcount.addFlops(2*A.nFlops());
      flopcount.addSiteFlops(2*Nc*Ns,s);
    }

    TF x; x[s]=zero;
 
    // now work out r= chi - Apsi = chi - r0
    r[s] = b; 

    Double r_sq = norm2(r,s);
    flopcount.addSiteFlops(4*Nc*Ns,s);


    QDPIO::cout << "Reliable CG: || r0 ||/|| b ||=" << sqrt(r_sq/chi_norm) << endl;


    Double rNorm = sqrt(r_sq);
    Double r0Norm = rNorm;
    Double maxrx = rNorm;
    Double maxrr = rNorm;
    bool updateR = false;
    bool updateX = false;
    

    // Now initialise v = p = 0
    TF p;
    Double a, c, d;

    // The iterations 
    for(k = 0; k < MaxCG && !convP; k++) { 
      if( k == 0 ) { 
	p[s] = r;
      }
      else { 
	Double beta = r_sq / c;
	RF br = beta;
	p[s] = r + br*p;  flopcount.addSiteFlops(4*Nc*Ns,s);
      }

      c = r_sq;

      TF mmp,mp;
      AF(mp, p, PLUS); 
      d = norm2(mp,s); 
      AF(mmp,mp,MINUS); 

      a = c/d;
      RF ar = a;
      x[s] += ar*p;  
      r[s] -= ar*mmp; 

      r_sq = norm2(r,s); 
      
      //      flopcount.addSiteFlops(4*Nc*Ns,s); <mp, mp>
      //      flopcount.addSiteFlops(4*Nc*Ns,s); x += a * p
      //      flopcount.addSiteFlops(4*Nc*Ns,s); r -= a * mm
      //      flopcount.addSiteFlops(4*Nc*Ns,s); norm2(r)
      flopcount.addSiteFlops(16*Nc*Ns,s);
      flopcount.addFlops(2*A.nFlops());

      // Reliable update part...
      rNorm = sqrt(r_sq);
      if( toBool( rNorm > maxrx) ) maxrx = rNorm;
      if( toBool( rNorm > maxrr) ) maxrr = rNorm;
      
      updateX = toBool ( rNorm < Delta*r0Norm && r0Norm <= maxrx );
      updateR = toBool ( rNorm < Delta*maxrr && r0Norm <= maxrr ) || updateX;

      // Do the R update with real DP residual
      if( updateR ) { 

	{
	  T tmp1,tmp2;
	  x_dble[s] = x;
	  
	  A(tmp1, x_dble, PLUS); // Use full solution so far
	  A(tmp2, tmp1, MINUS); // Use full solution so far

	  r_dble[s] = b - tmp2;
	}

	r[s] = r_dble;     // new R = b - Ax
	r_sq = norm2(r_dble,s);

	flopcount.addSiteFlops(6*Nc*Ns,s); // 4 from norm2, 2 from r=b-tmp2
	flopcount.addFlops(2*A.nFlops());

	rNorm = sqrt(r_sq);
	maxrr = rNorm;
	
	// Group wise x update
	if( updateX ) { 
	  if( ! updateR ) { x_dble[s]=x; } // if updateR then this is done already
	  psi[s] += x_dble; // Add on group accumulated solution in y
	  flopcount.addSiteFlops(2*Nc*Ns,s);

	  x[s] = zero; // zero y
	  b[s] = r_dble;
	  r0Norm = rNorm;
	  maxrx = rNorm;
	}

      }

      // Convergence check
      if( toBool(r_sq < rsd_sq ) ) {
	// We've converged.

	// if updateX true, then we have just updated psi
	// strictly x[s] should be zero, so it should be OK to add it
	// but why do the work if you don't need to
	x_dble[s] = x;
	psi[s]+=x_dble;
	flopcount.addSiteFlops(2*Nc*Ns,s);
	ret.resid = rNorm;
	ret.n_count = k;
	convP = true;
      }
      else { 
	convP = false;
      }

    }

    // Loop is finished. Report FLOP Count...
    swatch.stop();
    flopcount.report("reliable_invcg2", swatch.getTimeInSeconds());

    // Check for nonconvergence
    if( k >= MaxCG ) { 
      QDPIO::cout << "Nonconvergence: Reliable CG Failed to converge in " << MaxCG << " iterations " << endl;
      QDP_abort(1);
    }

    // Done
    END_CODE();
    return ret;
  }



SystemSolverResults_t
InvCGReliable(const LinearOperator<LatticeFermionF>& A,
	      const LatticeFermionF& chi,
	      LatticeFermionF& psi,
	      const Real& RsdCG, 
	      const Real& Delta,
	      int MaxCG)
  
{
  return RelInvCG_a<LatticeFermionF,LatticeFermionF, RealF>(A,A, chi, psi, RsdCG, Delta, MaxCG);
}

  // Pure double
SystemSolverResults_t
InvCGReliable(const LinearOperator<LatticeFermionD>& A,
		    const LatticeFermionD& chi,
		    LatticeFermionD& psi,
		    const Real& RsdCG, 
		    const Real& Delta,
	      int MaxCG)
{
  return RelInvCG_a<LatticeFermionD, LatticeFermionD, RealD>(A,A, chi, psi, RsdCG, Delta, MaxCG);
}

  // single double
SystemSolverResults_t
InvCGReliable(const LinearOperator<LatticeFermionD>& A,
		    const LinearOperator<LatticeFermionF>& AF,
		    const LatticeFermionD& chi,
		    LatticeFermionD& psi,
		    const Real& RsdCG, 
		    const Real& Delta,
	      int MaxCG)
{
  return RelInvCG_a<LatticeFermionD, LatticeFermionF, RealF>(A,AF, chi, psi, RsdCG, Delta, MaxCG);
}


}  // end namespace Chroma
