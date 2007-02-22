// $Id: inv_rel_cg1.cc,v 3.1 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/inv_rel_cg1.h"

namespace Chroma {

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 */
#undef CHROMA_INV_REL_CG1_RSD_CHK

template<typename T>
void InvRelCG1_a(const LinearOperator<T>& A,
		 const T& chi,
		 T& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count)
{
  START_CODE();

  const Subset& s = A.subset();

  Real chi_sq = Real(norm2(chi,s));
  Real rsd_sq = (RsdCG * RsdCG)*chi_sq;

  //
  //  r[0]  :=  Chi - A . Psi[0] 
    
  //                    
  //  r  :=  [ Chi  -  A . psi ]


  psi[s] = zero;

  T tmp1, tmp2;
  T r;

  r[s] = chi;
  Real inner_tol ;

  //  p[1]  :=  r[0]
  T p;
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
    END_CODE();
    return;
  }

  //
  //  FOR k FROM 1 TO MaxCG DO
  //
  LatticeFermion q;
  q=zero;

  for(int k = 1; k <= MaxCG; ++k) {

    // The sqrt(norm2(p)) part is taken care of by using relative residua
    // in the inner solve
    inner_tol = sqrt(rsd_sq)*sqrt(zeta);
 
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

#ifdef CHROMA_INV_REL_CG1_RSD_CHK
    {
      LatticeFermion rcheck;
      A(rcheck, psi, PLUS);
      rcheck[s] -= chi;
      QDPIO::cout << "InvCG1: inter " << k<< " || b - Ax ||^2 = " << norm2(rcheck,s) << " || r ||^2 = " << c << endl;
    }
#endif

    zeta += Double(1)/c;

    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    Real b = Real(c) / Real(cp);
    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */

    cp = c;	                /* 2 Nc Ns  flops */

    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      END_CODE();
      return;
    }
  }
  n_count = MaxCG;
  QDPIO::cerr << "Nonconvergence Warning: n_count =" << n_count << endl;
  END_CODE();
}


// Fix here for now
template<>
void InvRelCG1(const LinearOperator<LatticeFermion>& A,
	       const LatticeFermion& chi,
	       LatticeFermion& psi,
	       const Real& RsdCG, 
	       int MaxCG, 
	       int& n_count)
{
  InvRelCG1_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

}  // end namespace Chroma
