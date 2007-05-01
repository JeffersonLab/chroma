// -*- C++ -*-
// $Id: invmr.h,v 3.3 2007-05-01 15:27:27 bjoo Exp $
/*! \file
 *  \brief Minimal-Residual (MR) for a generic fermion Linear Operator
 */

#ifndef __invmr_h__
#define __invmr_h__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma 
{

  //! Minimal-residual (MR) algorithm for a generic Linear Operator 
  /*! \ingroup invert
   * This subroutine uses the Minimal Residual (MR) algorithm to determine
   * the solution of the set of linear equations. Here we allow M to be nonhermitian.
   *
   *   	    Chi  =  M . Psi 
   *
   * Algorithm:
   *
   *  Psi[0]                                      Argument
   *  r[0]    :=  Chi  -  M . Psi[0] ;            Initial residual
   *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;       Converged?
   *  FOR k FROM 1 TO MaxCG DO                    MR iterations
   *      a[k-1]  := <M.r[k-1],r[k-1]> / <M.r[k-1],M.r[k-1]> ;
   *      ap[k-1] := MRovpar * a[k] ;             Overrelaxtion step
   *      Psi[k]  += ap[k-1] r[k-1] ;             New solution vector
   *      r[k]    -= ap[k-1] A . r[k-1] ;         New residual
   *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;   Converged?

   * Arguments:

   *  \param M       Linear Operator             (Read)
   *  \param chi     Source                      (Read)
   *  \param psi     Solution                    (Modify)
   *  \param RsdCG   MR residual accuracy        (Read)
   *  \param MRovpar Overrelaxation parameter    (Read)
   *  \param MaxCG   Maximum MR iterations       (Read)

   * Local Variables:

   *  r   	Residual vector
   *  cp  	| r[k] |**2
   *  c   	| r[k-1] |**2
   *  k   	MR iteration counter
   *  a   	a[k]
   *  d   	< M.r[k], M.r[k] >
   *  R_Aux     Temporary for  M.Psi
   *  Mr        Temporary for  M.r

   * Global Variables:

   *  MaxCG       Maximum number of MR iterations allowed
   *  RsdCG       Maximum acceptable MR residual (relative to source)
   *
   * Subroutines:
   *
   *  M           Apply matrix to vector
   *
   * @{
   */

  template<typename T>
  SystemSolverResults_t 
  InvMR(const LinearOperator<T>& M,
	const T& chi,
	T& psi,
	const Real& MRovpar,
	const Real& RsdMR, 
	int MaxMR,
	enum PlusMinus isign);

  /*! @} */  // end of group invert

}  // end namespace Chroma

#endif
