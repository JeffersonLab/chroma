// $Id: invcg2_timing_hacks.cc,v 1.2 2004-03-22 15:29:29 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg2_timing_hacks.h"

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
void InvCG2_timing_hacks_a(const LinearOperator<T>& M,
			   const T& chi,
			   T& psi,
			   const Real& RsdCG, 
			   int MaxCG, 
			   int& n_count)
{
  const OrderedSubset& s = M.subset();

//  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));
  Real chi_sq =  Real(norm2(chi,s));

  // QDPIO::cout << "chi_norm = " << sqrt(chi_sq) << endl;
  Real rsd_sq = (RsdCG * RsdCG) * chi_sq;

  //                                            +
  //  r[0]  :=  Chi - A . Psi[0]    where  A = M  . M
    
  //                      +
  //  r  :=  [ Chi  -  M(u)  . M(u) . psi ]
  T  r, mp, mmp;
  M(mp, psi, PLUS);
  M(mmp, mp, MINUS);
  r[s] = chi - mmp;

  //  p[1]  :=  r[0]
  T p;
  p[s] = r;
  
  //  Cp = |r[0]|^2
  Double cp = norm2(r, s);   	       	   /* 2 Nc Ns  flops */

  // QDPIO::cout << "InvCG: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq << endl;
  
  // Disable early termination
#if 0
  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  if ( toBool(cp  <=  rsd_sq) )
  {
    n_count = 0;
    return;
  }
#endif

  //
  //  FOR k FROM 1 TO n_count  DO
  //
  Real a, b;
  Double c, d;
  
  // Go to n_count iters
  for(int k = 1; k <= n_count; ++k)
  {
    //  c  =  | r[k-1] |**2
    c = cp;

    //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
    //      	       	       	       	       	  +
    //  First compute  d  =  < p, A.p >  =  < p, M . M . p >  =  < M.p, M.p >
    //  Mp = M(u) * p
    M(mp, p, PLUS);

    //  d = | mp | ** 2
    d = norm2(mp, s);	/* 2 Nc Ns  flops */

    // Disable this. Set a = 1 to stop convergence
#if 0
    a = Real(c)/Real(d);
#else 
    a = Real(1)/Real(1);
#endif

    //  Psi[k] += a[k] p[k]
    psi[s] += a * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] A . p[k] ;
    //      	       +            +
    //  r  =  r  -  M(u)  . Mp  =  M  . M . p  =  A . p
    M(mmp, mp, MINUS);
    r[s] -= a * mmp;

    //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;

    //  cp  =  | r[k] |**2
    cp = norm2(r, s);	                /* 2 Nc Ns  flops */

    // QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << endl;

    // Disable termination
#if 0
    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      return;
    }
#endif

#if 0
    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    b = Real(cp) / Real(c);
#else
    b = Real(1)/Real(1);
#endif

    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */
  }
  // n_count = MaxCG;
  // QDP_error_exit("too many CG iterations: count = %d", n_count);
}


// Fix here for now
template<>
void InvCG2_timing_hacks(const LinearOperator<LatticeFermion>& M,
	    const LatticeFermion& chi,
	    LatticeFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count)
{
  InvCG2_timing_hacks_a(M, chi, psi, RsdCG, MaxCG, n_count);
}

// Fix here for now
template<>
void InvCG2_timing_hacks(const LinearOperator<LatticeDWFermion>& M,
	    const LatticeDWFermion& chi,
	    LatticeDWFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count)
{
  InvCG2_timing_hacks_a(M, chi, psi, RsdCG, MaxCG, n_count);
}
