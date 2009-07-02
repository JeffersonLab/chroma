// -*- C++ -*-
// $Id: syssolver_linop_bicrstab.h,v 3.1 2009-07-02 22:11:03 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by BICGSTAB
 */

#ifndef __syssolver_linop_bicrstab_h__
#define __syssolver_linop_bicrstab_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_bicgstab_params.h"
#include "actions/ferm/invert/invbicrstab.h"

namespace Chroma
{

  //! BICGSTAB system solver namespace
  namespace LinOpSysSolverBiCRStabEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by BICGSTAB
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverBiCRStab : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverBiCRStab(Handle< LinearOperator<T> > A_,
		     const SysSolverBiCGStabParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverBiCRStab() {}

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
      res = InvBiCRStab(*A, 
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
    LinOpSysSolverBiCRStab() {}

    Handle< LinearOperator<T> > A;
    SysSolverBiCGStabParams invParam;
  };

} // End namespace

#endif 

