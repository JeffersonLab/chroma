// $Id: invcg2_timing_hacks.cc,v 3.2 2007-02-22 21:11:46 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg2_timing_hacks.h"

namespace Chroma {

// This is a hack version of invcg2 designed to allow timing.
// In particular stuff that would make the CG converge has been
// DISABLED... The length of the CG is completely controlled by 
// the INPUT n_count

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 *
 * where       A = M^dag . M
 *
 * Algorithm:

 *  Psi[0]  :=  initial guess;    	       Linear interpolation (argument)
 *  r[0]    :=  Chi - M^dag . M . Psi[0] ;     Initial residual
 *  p[1]    :=  r[0] ;	       	       	       Initial direction
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;      Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       CG iterations
 *      a[k] := |r[k-1]|**2 / <Mp[k],Mp[k]> ;
 *      Psi[k] += a[k] p[k] ;   	       New solution vector
 *      r[k] -= a[k] M^dag . M . p[k] ;        New residual
 *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;  Converged?
 *      b[k+1] := |r[k]|**2 / |r[k-1]|**2 ;
 *      p[k+1] := r[k] + b[k+1] p[k];          New direction
 *
 * Arguments:
 *
 *  \param M       Linear Operator    	       (Read)
 *  \param chi     Source	               (Read)
 *  \param psi     Solution    	    	       (Modify)
 *  \param RsdCG   CG residual accuracy        (Read)
 *  \param MaxCG   Maximum CG iterations       (Read)
 *  \param n_count Number of CG iteration      (Write)
 *
 * Local Variables:
 *
 *  p   	       Direction vector
 *  r   	       Residual vector
 *  cp  	       | r[k] |**2
 *  c   	       | r[k-1] |**2
 *  k   	       CG iteration counter
 *  a   	       a[k]
 *  b   	       b[k+1]
 *  d   	       < p[k], A.p[k] >
 *  Mp  	       Temporary for  M.p
 *
 * Subroutines:
 *                             +               
 *  A       Apply matrix M or M  to vector
 *
 * Operations:
 *
 *  2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns )
 */

template<typename T>
SystemSolverResults_t
InvCG2_timings_a(const LinearOperator<T>& M,
		 const T& chi,
		 T& psi,
		 int MaxCG)
			
{
    START_CODE();

    const Subset& s = M.subset();

    SystemSolverResults_t  res;
    T mp;                moveToFastMemoryHint(mp);
    T mmp;               moveToFastMemoryHint(mmp);
    T p;                 moveToFastMemoryHint(p);
    moveToFastMemoryHint(psi,true);
    T r;                 moveToFastMemoryHint(r);
    T chi_internal;      moveToFastMemoryHint(chi_internal);

    chi_internal[s] = chi;

    QDPIO::cout << "InvCG2: starting" << endl;
    FlopCounter flopcount;
    flopcount.reset();
    StopWatch swatch;
    swatch.reset();
    swatch.start();


    Real chi_sq =  Real(norm2(chi_internal,s));
    flopcount.addSiteFlops(4*Nc*Ns,s);


    //                                            +
    //  r[0]  :=  Chi - A . Psi[0]    where  A = M  . M
    
    //                      +
    //  r  :=  [ Chi  -  M(u)  . M(u) . psi ]
    M(mp, psi, PLUS);
    M(mmp, mp, MINUS);
    flopcount.addFlops(2*M.nFlops());

    r[s] = chi_internal - mmp;
    flopcount.addSiteFlops(2*Nc*Ns,s);

    //  p[1]  :=  r[0]
    p[s] = r;
  
    //  Cp = |r[0]|^2
    Double cp = norm2(r, s);   	       	   /* 2 Nc Ns  flops */
    flopcount.addSiteFlops(4*Nc*Ns, s);

    //
    //  FOR k FROM 1 TO MaxCG DO
    //
    Real a, b;
    Double c, d;
    int k;
    for(k = 1; k <= MaxCG; ++k)
    {
      //  c  =  | r[k-1] |**2
      c = cp;

      //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
      //      	       	       	       	       	  +
      //  First compute  d  =  < p, A.p >  =  < p, M . M . p >  =  < M.p, M.p >
      //  Mp = M(u) * p
      M(mp, p, PLUS);  flopcount.addFlops(M.nFlops());

      //  d = | mp | ** 2
      d = norm2(mp, s);  flopcount.addSiteFlops(4*Nc*Ns,s);

      a = Real(1)/Real(1);

      //  Psi[k] += a[k] p[k]
      psi[s] += a * p;    flopcount.addSiteFlops(4*Nc*Ns,s);

      //  r[k] -= a[k] A . p[k] ;
      //      	       +            +
      //  r  =  r  -  M(u)  . Mp  =  M  . M . p  =  A . p
      M(mmp, mp, MINUS);
      flopcount.addFlops(M.nFlops());

      r[s] -= a * mmp;
      flopcount.addSiteFlops(4*Nc*Ns, s);

      //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;

      //  cp  =  | r[k] |**2
      cp = norm2(r, s);    flopcount.addSiteFlops(4*Nc*Ns,s);

      // Disable convergence
      b = Real(1)/Real(1);

      //  p[k+1] := r[k] + b[k+1] p[k]
      p[s] = r + b*p;    flopcount.addSiteFlops(4*Nc*Ns,s);
    }
    
    
    swatch.stop();
    flopcount.report("invcg2_timings", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);
    res.n_count = MaxCG;
    res.resid   = 1.0e-30;
    QDPIO::cout << "InvCG2 Timings: Iteration Count = " << k << endl;
    END_CODE();
    return res;


}


// Fix here for now
template<>
SystemSolverResults_t
InvCG2_timings(const LinearOperator<LatticeFermion>& M,
	       const LatticeFermion& chi,
	       LatticeFermion& psi,
	       int MaxCG)
{
  return InvCG2_timings_a(M, chi, psi, MaxCG);
}

}  // end namespace Chroma
