// -*- C++ -*-
// $Id: syssolver_linop_bicgstab.h,v 3.2 2008-04-05 19:04:38 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by BICGSTAB
 */

#ifndef __syssolver_linop_bicgstab_h__
#define __syssolver_linop_bicgstab_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_bicgstab_params.h"
#include "actions/ferm/invert/invbicgstab.h"

namespace Chroma
{

  //! BICGSTAB system solver namespace
  namespace LinOpSysSolverBiCGStabEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by BICGSTAB
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverBiCGStab : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverBiCGStab(Handle< LinearOperator<T> > A_,
		     const SysSolverBiCGStabParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverBiCGStab() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
    {
      START_CODE();
      
      
      SystemSolverResults_t res;  // initialized by a constructor
      
      // For now solve with PLUS until we add a way to explicitly
      // ask for MINUS
      res = InvBiCGStab(*A, 
			chi, 
			psi, 
			invParam.RsdBiCGStab, 
			invParam.MaxBiCGStab, 
			PLUS);
      
      
      END_CODE();
      
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverBiCGStab() {}

    Handle< LinearOperator<T> > A;
    SysSolverBiCGStabParams invParam;
  };

} // End namespace

#endif 

