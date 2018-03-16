// -*- C++ -*-
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invcg_adj__
#define __invcg_adj__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma 
{

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
   *      Psi[k] += a[k] p[k] ;   	       New solution std::vector
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
   *  p   	       Direction std::vector
   *  r   	       Residual std::vector
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
   *  A       Apply matrix M or M  to std::vector
   *
   * Operations:
   *
   *  2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns )
   *
   * @{
   */

  Double norm2_adj(const LatticeColorMatrix& C,const Subset& s){
    return Double(0.5)*real(sum(trace(adj(C)*C),s)) ;
  }

  DComplex inner_prod_adj(const LatticeColorMatrix& C,
			const LatticeColorMatrix& X,
			const Subset& s){
    return Double(0.5)*sum(trace(adj(C)*X),s) ;
  }

  Double norm2_adj(const LatticeColorMatrix& C){
    return Double(0.5)*real(sum(trace(adj(C)*C))) ;
  }
  
  DComplex inner_prod_adj(const LatticeColorMatrix& C,
			  const LatticeColorMatrix& X){
    return Double(0.5)*sum(trace(adj(C)*X)) ;
  }
  
  // Single precision
  SystemSolverResults_t 
  InvCG_adj(const LinearOperator<LatticeColorMatrix>& M,
	 const LatticeColorMatrix& chi,
	 LatticeColorMatrix& psi,
	 const Real& RsdCG, 
	 int MaxCG);


  /*! @} */  // end of group invert

}  // end namespace Chroma

#endif
