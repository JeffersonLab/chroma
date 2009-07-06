// -*- C++ -*-
// $Id: reliable_ibicgstab.h,v 3.1 2009-07-06 19:02:34 bjoo Exp $
/*! \file
 *  \brief BiCGStab Solver with reliable updates
 */

#ifndef __reliable_invibicgstab__
#define __reliable_invibicgstab__

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
  InvIBiCGStabReliable(const LinearOperator<LatticeFermionF>& A,
		      const LatticeFermionF& chi,
		      LatticeFermionF& psi,
		      const Real& RsdBiCGStab, 
		      const Real& Delta,
		      int MaxBiCGStab, 
		      enum PlusMinus isign);

  // Pure double
  SystemSolverResults_t
  InvIBiCGStabReliable(const LinearOperator<LatticeFermionD>& A,
		      const LatticeFermionD& chi,
		      LatticeFermionD& psi,
		      const Real& RsdBiCGStab, 
		      const Real& Delta,
		      int MaxBiCGStab, 
		      enum PlusMinus isign);

   // single double
  SystemSolverResults_t
  InvIBiCGStabReliable(const LinearOperator<LatticeFermionD>& A,
		      const LinearOperator<LatticeFermionF>& AF,
		      const LatticeFermionD& chi,
		      LatticeFermionD& psi,
		    const Real& RsdBiCGStab, 
		      const Real& Delta,
		      int MaxBiCGStab, 
		      enum PlusMinus isign);
 

  /*! @} */  // end of group invert
	    
}  // end namespace Chroma

#endif
