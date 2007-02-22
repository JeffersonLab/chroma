// $Id: inv_rel_cg2.cc,v 3.1 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/inv_rel_cg2.h"

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
void InvRelCG2_a(const LinearOperator<T>& M,
		 const T& chi,
		 T& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count)
{
  const Subset& s = M.subset();
  
  //  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));
  Real chi_sq =  Real(norm2(chi,s));

  QDPIO::cout << "chi_norm = " << sqrt(chi_sq) << endl;
  Real rsd_sq = (RsdCG * RsdCG) * chi_sq;

  //                                            +
  //  r[0]  :=  Chi - A . Psi[0]    where  A = M  . M
    
  //                      +
  //  r  :=  [ Chi  -  M(u)  . M(u) . psi ]
  T  r, mp, mmp;
  psi[s] = zero;

  r[s] = chi;
  //  p[1]  :=  r[0]
  T p;
  p = zero;
  p[s] = r;
  
  //  Cp = |r[0]|^2
  Double c = norm2(r, s);   	       	   /* 2 Nc Ns  flops */
  Double cp = c;

  Double zeta = Double(1)/c;

  QDPIO::cout << "InvRelCG2: k = 0  c = " << cp << "  rsd_sq = " << rsd_sq << endl;

  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  if ( toBool(c  <=  rsd_sq) )
  {
    n_count = 0;
    return;
  }

  //
  //  FOR k FROM 1 TO MaxCG DO
  //
  Real a;
  Double d;
  
  for(int k = 1; k <= MaxCG; ++k)
  {
    // Inner tolerance = epsilon || chi || || p || sqrt(zeta) / 2
    //
    // The || p || part is taken care of the fact that we are using 
    // relative residua in the inner solve. The factor of 2 is because
    // we apply the operator twice...
    Real inner_tol = sqrt(rsd_sq)*sqrt(zeta)/Real(2);

    // Compute M^{dag} M p
    M(mp, p, PLUS, inner_tol);
    M(mmp, mp, MINUS, inner_tol);

    //  d = < M^{dag} M p, p>
    d = innerProductReal(mmp, p, s);	/* 2 Nc Ns  flops */

    //  a[k] := | r |**2 / < p[k], Ap[k] > ;
    a = Real(c)/Real(d);

    //  Psi[k] += a[k] p[k]
    psi[s] += a * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] M^{dag}M. p[k] ;
    r[s] -= a * mmp;

    //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;

    //  cp  =  | r[k] |**2
    c = norm2(r, s);	                /* 2 Nc Ns  flops */

    // update relaxation factor
    zeta += Double(1)/c;

    // Update p as:
    // p[k+1] := r[k] + c/cp * p
    //
    // we put c/cp into b
    //
    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    Real b = Real(c) / Real(cp);

    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */

    cp = c;

    QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << endl;

    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      return;
    }
  }
  n_count = MaxCG;
  QDPIO::cerr << "Nonconvergence Warning n_count = " << n_count << endl;
}


// Fix here for now
template<>
void InvRelCG2(const LinearOperator<LatticeFermion>& M,
	       const LatticeFermion& chi,
	       LatticeFermion& psi,
	       const Real& RsdCG, 
	       int MaxCG, 
	       int& n_count)
{
  InvRelCG2_a(M, chi, psi, RsdCG, MaxCG, n_count);
}

}  // end namespace Chroma
