// -*- C++ -*-
// $Id: reliable_cg.h,v 3.1 2009-05-22 14:21:39 bjoo Exp $
/*! \file
 *  \brief BiCGStab Solver with reliable updates
 */

#ifndef __reliable_cg_h__
#define __reliable_cg_h__

#include "linearop.h"
#include "syssolver.h"

namespace Chroma 
{

  //! Bi-CG stabilized
  /*! \ingroup invert
   *
   * @{
   */
  SystemSolverResults_t
  InvCGReliable(const LinearOperator<LatticeFermionF>& A,
		const LatticeFermionF& chi,
		LatticeFermionF& psi,
		const Real& RsdCG, 
		const Real& Delta,
		int MaxCG);

  // Pure double
  SystemSolverResults_t
  InvCGReliable(const LinearOperator<LatticeFermionD>& A,
		const LatticeFermionD& chi,
		LatticeFermionD& psi,
		const Real& RsdCG, 
		const Real& Delta,
		int MaxCG);

  
   // single double
  SystemSolverResults_t
  InvCGReliable(const LinearOperator<LatticeFermionD>& A,
		const LinearOperator<LatticeFermionF>& AF,
		const LatticeFermionD& chi,
		LatticeFermionD& psi,
		const Real& RsdCG, 
		const Real& Delta,
		int MaxCG);
  

  /*! @} */  // end of group invert
	    
}  // end namespace Chroma

#endif
