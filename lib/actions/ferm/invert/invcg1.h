// -*- C++ -*-
// $Id: invcg1.h,v 3.2 2008-03-25 10:46:51 mcneile Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invcg1__
#define __invcg1__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma {

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 *
 * where       A is hermitian
 *
 * Algorithm:

 *  Psi[0]  :=  initial guess;    	       Linear interpolation (argument)
 *  r[0]    :=  Chi - A. Psi[0] ;     Initial residual
 *  p[1]    :=  r[0] ;	       	       	       Initial direction
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;      Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       CG iterations
 *      a[k] := |r[k-1]|**2 / <p[k],A p[k]> ;
 *      Psi[k] += a[k] p[k] ;   	       New solution vector
 *      r[k] -= a[k] A p[k] ;        New residual
 *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;  Converged?
 *      b[k+1] := |r[k]|**2 / |r[k-1]|**2 ;
 *      p[k+1] := r[k] + b[k+1] p[k];          New direction
 *
 * Arguments:
 *
 *  \param A       Linear Operator    	       (Read)
 *  \param chi     Source	               (Read)
 *  \param psi     Solution    	    	       (Modify)
 *  \param RsdCG   CG residual accuracy        (Read)
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
 *  ap  	       Temporary for  M.p
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
InvCG1(const LinearOperator<T>& A,
       const T& chi,
       T& psi,
       const Real& RsdCG, 
       int MaxCG, int MinCG=0);

}  // end namespace Chroma

#endif
