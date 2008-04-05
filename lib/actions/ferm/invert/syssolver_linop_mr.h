// -*- C++ -*-
// $Id: syssolver_linop_mr.h,v 1.4 2008-04-05 19:04:38 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by MR
 */

#ifndef __syssolver_linop_mr_h__
#define __syssolver_linop_mr_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_mr_params.h"
#include "actions/ferm/invert/invmr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace LinOpSysSolverMREnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by MR
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverMR : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverMR(Handle< LinearOperator<T> > A_,
		     const SysSolverMRParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverMR() {}

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
	{
	  res = InvMR(*A, chi, psi, invParam.MROver, invParam.RsdMR, invParam.MaxMR, PLUS);
	}

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverMR() {}

    Handle< LinearOperator<T> > A;
    SysSolverMRParams invParam;
  };

} // End namespace

#endif 

