// -*- C++ -*-
// $Id: invcg2_array.h,v 3.3 2007-08-27 18:18:54 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invcg2_array__
#define __invcg2_array__

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
   *
   * @{
   */
  // Single precision
  SystemSolverResults_t 
  InvCG2(const LinearOperatorArray<LatticeFermionF>& M,
	 const multi1d<LatticeFermionF>& chi,
	 multi1d<LatticeFermionF>& psi,
	 const Real& RsdCG, 
	 int MaxCG);

  // Double precision
  SystemSolverResults_t 
  InvCG2(const LinearOperatorArray<LatticeFermionD>& M,
	 const multi1d<LatticeFermionD>& chi,
	 multi1d<LatticeFermionD>& psi,
	 const Real& RsdCG, 
	 int MaxCG);

  /*! @} */  // end of group invert

}  // end namespace Chroma

#endif
