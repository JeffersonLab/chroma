// $Id: invcg2.cc,v 1.13 2005-01-14 20:13:05 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg2.h"

namespace Chroma {

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
void InvCG2_a(const LinearOperator<T>& M,
	      const T& chi,
	      T& psi,
	      const Real& RsdCG, 
	      int MaxCG, 
	      int& n_count)
{
  START_CODE();

  const OrderedSubset& s = M.subset();

//  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));
  Real chi_sq =  Real(norm2(chi,s));

#if 0
  QDPIO::cout << "chi_norm = " << sqrt(chi_sq) << endl;
#endif

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

#if 0
  QDPIO::cout << "InvCG: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq << endl;
#endif

  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  if ( toBool(cp  <=  rsd_sq) )
  {
    n_count = 0;
    END_CODE();
    return;
  }

  //
  //  FOR k FROM 1 TO MaxCG DO
  //
  Real a, b;
  Double c, d;
  
  for(int k = 1; k <= MaxCG; ++k)
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

    a = Real(c)/Real(d);

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

#if 0
    QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << endl;
#endif

    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      END_CODE();
      return;
    }

    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    b = Real(cp) / Real(c);

    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */
  }
  n_count = MaxCG;
  QDP_error_exit("too many CG iterations: count = %d", n_count);
}


// Fix here for now
template<>
void InvCG2(const LinearOperator<LatticeFermion>& M,
	    const LatticeFermion& chi,
	    LatticeFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count)
{
  InvCG2_a(M, chi, psi, RsdCG, MaxCG, n_count);
}


template<>
void InvCG2(const LinearOperator<LatticeStaggeredFermion>& M,
	    const LatticeStaggeredFermion& chi,
	    LatticeStaggeredFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count)
{
  InvCG2_a(M, chi, psi, RsdCG, MaxCG, n_count);
}

}  // end namespace Chroma
