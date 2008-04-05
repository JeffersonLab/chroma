// -*- C++ -*-
// $Id: syssolver_mdagm_mr.h,v 1.2 2008-04-05 19:04:38 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by MR
 */

#ifndef __syssolver_mdagm_mr_h__
#define __syssolver_mdagm_mr_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "ldag.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/syssolver_mr_params.h"
#include "actions/ferm/invert/invmr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace MdagMSysSolverMREnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a MR system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverMR : public MdagMSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverMR(Handle< LinearOperator<T> > A_,
		     const SysSolverMRParams& invParam_) : 
      A(A_), Adag(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMSysSolverMR() {}

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

	T tmp;      moveToFastMemoryHint(tmp);

	res = InvMR(*Adag, chi, tmp, invParam.RsdMR, invParam.MaxMR);
	res = InvMR(*A, tmp, psi, invParam.RsdMR, invParam.MaxMR);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMSysSolverMR() {}

    Handle< LinearOperator<T> > A;
    Handle< LinearOperator<T> > Adag;
    SysSolverMRParams invParam;
  };


} // End namespace

#endif 

