// $Id: inv_rel_cg1.cc,v 1.3 2004-05-18 12:40:15 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/inv_rel_cg1.h"

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 */

template<typename T>
void InvRelCG1_a(const ApproxLinearOperator<T>& A,
		 const T& chi,
		 T& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count)
{
  const OrderedSubset& s = A.subset();

  Real chi_sq = Real(norm2(chi,s));
  Real rsd_sq = (RsdCG * RsdCG)*chi_sq;

  //
  //  r[0]  :=  Chi - A . Psi[0] 
    
  //                    
  //  r  :=  [ Chi  -  A . psi ]



  T tmp1, tmp2;
  tmp1 = tmp2 = zero; 
  T r = zero;
  Real inner_tol ;

  // A is Hermitian
  A(tmp1, psi, PLUS);
  r[s] = chi - tmp1;

  //  p[1]  :=  r[0]
  T p = zero;
  p[s] = r;

  
  //  Cp = |r[0]|^2
  Double c = norm2(r, s);   	       	   /* 2 Nc Ns  flops */
  Double cp = c;
  Double zeta = Double(1)/c;

  QDPIO::cout << "InvRelCG1: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq 
<< endl;

  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  if ( toBool(cp  <=  rsd_sq) )
  {
    n_count = 0;
    return;
  }

  //
  //  FOR k FROM 1 TO MaxCG DO
  //
  LatticeFermion q;
  q=zero;

  for(int k = 1; k <= MaxCG; ++k) {

    // Why do I need the fudge factor 100?
    inner_tol = sqrt(rsd_sq)*sqrt(norm2(p))*sqrt(zeta);
 
    A(q,p,PLUS,inner_tol);
  
#if 0
    // Funkyness work out || Ap - q || / || p || 
    LatticeFermion tmp;
    tmp = zero;
    A(tmp,p,PLUS);
    tmp[s] -= q;
    Double tmpnorm=sqrt(norm2(tmp,s));
    QDPIO::cout << "Inner tol = " << inner_tol << " || Ap - q ||  = " << tmpnorm << endl;
#endif 

    //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
    //      	       	       	       	       	  +
    //  First compute  d  =  < q,p > 
    Double d = innerProductReal(q, p, s);   
    //    a = Real(c);
    Real a = Real(c)/Real(d);

    //  Psi[k] += a[k] p[k]
    psi[s] += a * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] A . p[k] ;
    //      	      
    //  r  =  r  - a A p  
    r[s] -= a * q;

    c = norm2(r,s);
    zeta += Double(1)/c;

    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    Real b = Real(c) / Real(cp);
    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */

    cp = c;	                /* 2 Nc Ns  flops */

    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      return;
    }
  }
  n_count = MaxCG;
  QDP_error_exit("too many CG iterations: count = %d", n_count);
}


// Fix here for now
template<>
void InvRelCG1(const ApproxLinearOperator<LatticeFermion>& A,
	       const LatticeFermion& chi,
	       LatticeFermion& psi,
	       const Real& RsdCG, 
	       int MaxCG, 
	       int& n_count)
{
  InvRelCG1_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

