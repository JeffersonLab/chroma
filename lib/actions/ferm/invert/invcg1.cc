// $Id: invcg1.cc,v 3.5 2008-03-25 10:46:51 mcneile Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg1.h"

using namespace QDP::Hints;

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
 *  \return res    System solver results
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
SystemSolverResults_t 
InvCG1_a(const LinearOperator<T>& A,
	 const T& chi,
	 T& psi,
	 const Real& RsdCG, 
	 int MaxCG, int MinCG=0)
{
  START_CODE();

  const Subset& s = A.subset();

  SystemSolverResults_t  res;

  T p;                 moveToFastMemoryHint(p);
  T ap;                moveToFastMemoryHint(ap);
  
  // Hint for psi to be moved with copy
  moveToFastMemoryHint(psi,true);

  T r;                 moveToFastMemoryHint(r);
  T chi_internal;      moveToFastMemoryHint(chi_internal);

  // Can't move chi to fast space, because its const.
  chi_internal[s]=chi;

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
    res.n_count = 0;
    res.resid   = sqrt(cp);
    revertFromFastMemoryHint(psi, true); // Revert psi, and copy contents
    END_CODE();
    return res;
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

    if ( toBool(cp  <=  rsd_sq) && (toBool(MinCG <= 0 ) || toBool( k >= MinCG  )))
    {
      res.n_count = k;
      res.resid   = sqrt(cp);
      
      // Compute the actual residual
      {
	A(ap, psi, PLUS);
	Double actual_res = norm2(chi - ap,s);
	res.resid = sqrt(actual_res);
#if 0
    QDPIO::cout << "InvCG1:  true residual:  " << res.resid << endl;
#endif
      }
      revertFromFastMemoryHint(psi,true); // Get back psi, copy contents.
      END_CODE();
      return res;
    }

    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    b = Real(cp) / Real(c);

    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + b*p;	/* Nc Ns  flops */
  }
  res.n_count = MaxCG;
  res.resid   = sqrt(cp);
  QDP_error_exit("too many CG iterations: count = %d", res.n_count);
  revertFromFastMemoryHint(psi,true);
  END_CODE();
  return res;
}

// Fix here for now
template<>
SystemSolverResults_t 
InvCG1(const LinearOperator<LatticeFermion>& A,
       const LatticeFermion& chi,
       LatticeFermion& psi,
       const Real& RsdCG, 
       int MaxCG,int MinCG)
{
  return InvCG1_a(A, chi, psi, RsdCG, MaxCG,MinCG);
}

template<>
SystemSolverResults_t 
InvCG1(const LinearOperator<LatticeStaggeredFermion>& A,
       const LatticeStaggeredFermion& chi,
       LatticeStaggeredFermion& psi,
       const Real& RsdCG, 
       int MaxCG,int MinCG)
{
  return InvCG1_a(A, chi, psi, RsdCG, MaxCG,MinCG);
}

}  // end namespace Chroma
