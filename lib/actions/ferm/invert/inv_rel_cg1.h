// -*- C++ -*-
// $Id: inv_rel_cg1.h,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __inv_rel_cg1__
#define __inv_rel_cg1__

#include "linearop.h"

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
 *  c := cp := || r[0] ||^2
 *  zeta    := 1/c;
 *
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;      Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       CG iterations
 *      inner_tol := RsdCG*||chi||*||p||*sqrt(zeta)
 *
 *      q    := A(inner_tol) p
 *
 *      a[k] := |r[k-1]|**2 / <q,p[k]> ;
 *      Psi[k] += a[k] p[k] ;   	       New solution vector
 *      r[k] -= a[k] q ;        New residual
 *      c := <r , r>
 *      zeta := zeta + 1/c
 *
 *      b[k+1] := |r[k]|**2 / |r[k-1]|**2 = c/cp;
 *      p[k+1] := r[k] + b[k+1] p[k];          New direction
 *
 *      cp := c
 *      IF ( c < RsdCG^2 || chi || ) RETURN
 *
 * Arguments:
 *
 *  \param A       Linear Operator    	       (Read)
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
void InvRelCG1(const LinearOperator<T>& A,
	       const T& chi,
	       T& psi,
	       const Real& RsdCG, 
	       int MaxCG, 
	       int& n_count);

}  // end namespace Chroma


#endif
