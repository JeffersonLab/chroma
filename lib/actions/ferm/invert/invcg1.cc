// $Id: invcg1.cc,v 1.12 2005-06-27 18:06:32 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg1.h"

namespace Chroma {

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 *
 * where       A  is Hermitian Positive Definite
 *
 * Algorithm:

 *  Psi[0]  :=  initial guess;    	       Linear interpolation (argument)
 *  r[0]    :=  Chi - A. Psi[0] ;              Initial residual
 *  p[1]    :=  r[0] ;	       	       	       Initial direction
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;      Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       CG iterations
 *      a[k] := |r[k-1]|**2 / <p[k],Ap[k]> ;
 *      Psi[k] += a[k] p[k] ;   	       New solution vector
 *      r[k] -= a[k] A. p[k] ;        New residual
 *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;  Converged?
 *      b[k+1] := |r[k]|**2 / |r[k-1]|**2 ;
 *      p[k+1] := r[k] + b[k+1] p[k];          New direction
 *
 * Arguments:
 *
 *  \param M       Linear Operator    	       (Read)
 *  \param chi     Source	               (Read)
 *  \param psi     Solution    	    	       (Modify)
 *  \param RsdCG   CG residual accuracy        (Rea/Write)
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
 *  ap  	       Temporary for  A.p
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
void InvCG1_a(const LinearOperator<T>& A,
	      const T& chi,
	      T& psi,
	      const Real& RsdCG, 
	      int MaxCG, 
	      int& n_count)
{
  START_CODE();

  const OrderedSubset& s = A.subset();

  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));

  //
  //  r[0]  :=  Chi - A . Psi[0] 
    
  //                    
  //  r  :=  [ Chi  -  A . psi ]

  T tmp1, tmp2;
  tmp1 = tmp2 = zero; 
  T r = zero;

  // A is Hermitian
  A(tmp1, psi, PLUS);


  r[s] = chi - tmp1;

  //  p[1]  :=  r[0]
  T p = zero;
  p[s] = r;
  
  //  Cp = |r[0]|^2
  Double cp = norm2(r, s);   	       	   /* 2 Nc Ns  flops */

#if 0
  QDPIO::cout << "InvCG1: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq 
<< endl;
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
  T    ap;
  Real b;
  Double c;
  Complex a;
  Real d;

  for(int k = 1; k <= MaxCG; ++k)
  {
    //  c  =  | r[k-1] |**2
    c = cp;

    //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
    //      	       	       	       	       	  +
    //  First compute  d  =  < p, A.p > 
   

    // SCRAP THIS IDEA FOR NOW AND DO <p, A.p> TO KEEP TRACK OF 
    // "PRECONDITIONING" So ap =MdagMp

    A(ap, p, PLUS);
    
    // d = <p,A.p>
    d = innerProductReal(p, ap, s);

    //    a = Real(c);
    a = Real(c)/d;


    //  Psi[k] += a[k] p[k]
    psi[s] += a * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] A . p[k] ;
    //      	      
    //  r  =  r  - a A p  

    r[s] -= a * ap;



    //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;

    //  cp  =  | r[k] |**2
    cp = norm2(r, s);	                /* 2 Nc Ns  flops */

#if 0
    QDPIO::cout << "InvCG1: k = " << k << "  cp = " << cp << endl;
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


void InvCG1_a(const LinearOperator<LatticeFermion>& A,
	      const LatticeFermion& chi,
	      LatticeFermion& psi,
	      const Real& RsdCG, 
	      int MaxCG, 
	      int& n_count)
{
  START_CODE();

  const OrderedSubset& s = A.subset();

  LatticeFermion ap;                     ap.moveToFastMemoryHint();
  LatticeFermion r;                      r.moveToFastMemoryHint();
  LatticeFermion p;                      p.moveToFastMemoryHint();
  LatticeFermion chi_internal;           chi_internal.moveToFastMemoryHint();

  chi_internal[s] = chi;
 
  psi.moveToFastMemoryHint(true);


  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi_internal,s));

  //
  //  r[0]  :=  Chi - A . Psi[0] 
    
  //                    
  //  r  :=  [ Chi  -  A . psi ]

  // A is Hermitian
  A(ap, psi, PLUS);


  r[s] = chi_internal - ap;

  //  p[1]  :=  r[0]
  p[s] = r;
  
  //  Cp = |r[0]|^2
  Double cp = norm2(r, s);   	       	   /* 2 Nc Ns  flops */

#if 0
  QDPIO::cout << "InvCG1: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq 
<< endl;
#endif

  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  if ( toBool(cp  <=  rsd_sq) )
  {
    n_count = 0;
    psi.revertFromFastMemoryHint(true);
    END_CODE();
    return;
  }

  //
  //  FOR k FROM 1 TO MaxCG DO
  //
  
  Real b;
  Double c;
  Complex a;
  Real d;

  for(int k = 1; k <= MaxCG; ++k)
  {
    //  c  =  | r[k-1] |**2
    c = cp;

    //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
    //      	       	       	       	       	  +
    //  First compute  d  =  < p, A.p > 
   

    // SCRAP THIS IDEA FOR NOW AND DO <p, A.p> TO KEEP TRACK OF 
    // "PRECONDITIONING" So ap =MdagMp

    A(ap, p, PLUS);
    
    // d = <p,A.p>
    d = innerProductReal(p, ap, s);

    //    a = Real(c);
    a = Real(c)/d;


    //  Psi[k] += a[k] p[k]
    psi[s] += a * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] A . p[k] ;
    //      	      
    //  r  =  r  - a A p  

    r[s] -= a * ap;



    //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;

    //  cp  =  | r[k] |**2
    cp = norm2(r, s);	                /* 2 Nc Ns  flops */

#if 0
    QDPIO::cout << "InvCG1: k = " << k << "  cp = " << cp << endl;
#endif

    if ( toBool(cp  <=  rsd_sq) )
    {
      n_count = k;
      psi.revertFromFastMemoryHint(true);
      END_CODE();
      return;
    }

    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    b = Real(cp) / Real(c);

    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */
  }
  n_count = MaxCG;
  QDPIO::cerr << "nonconvergence warning: iters = " << n_count << endl;
  psi.revertFromFastMemoryHint(true);
  END_CODE();

}


// Fix here for now
template<>
void InvCG1(const LinearOperator<LatticeFermion>& A,
	    const LatticeFermion& chi,
	    LatticeFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count)
{
  InvCG1_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

template<>
void InvCG1(const LinearOperator<LatticeStaggeredFermion>& A,
	    const LatticeStaggeredFermion& chi,
	    LatticeStaggeredFermion& psi,
	    const Real& RsdCG, 
	    int MaxCG, 
	    int& n_count)
{
  InvCG1_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

}  // end namespace Chroma
