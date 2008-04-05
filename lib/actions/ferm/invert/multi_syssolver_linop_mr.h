// -*- C++ -*-
// $Id: multi_syssolver_linop_mr.h,v 1.2 2008-04-05 19:04:38 edwards Exp $
/*! \file
 *  \brief Solve a (M+shift)*psi=chi linear system by MR
 */

#ifndef __multi_syssolver_linop_mr_h__
#define __multi_syssolver_linop_mr_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_linop.h"
#include "actions/ferm/invert/multi_syssolver_mr_params.h"
#include "actions/ferm/invert/minvmr.h"


namespace Chroma
{

  //! MR system solver namespace
  namespace LinOpMultiSysSolverMREnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a MR system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpMultiSysSolverMR : public LinOpMultiSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpMultiSysSolverMR(Handle< LinearOperator<T> > A_,
			  const MultiSysSolverMRParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpMultiSysSolverMR() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<Real>& shifts, const T& chi) const
      {
	START_CODE();

	multi1d<Real> RsdCG(shifts.size());
	if (invParam.RsdCG.size() == 1)
	{
	  RsdCG = invParam.RsdCG[0];
	}
	else if (invParam.RsdCG.size() == RsdCG.size())
	{
	  RsdCG = invParam.RsdCG;
	}
	else
	{
	  QDPIO::cerr << "LinOpMultiSysSolverMR: shifts incompatible" << endl;
	  QDP_abort(1);
	}

	SystemSolverResults_t res;
  	MInvMR(*A, chi, psi, shifts, RsdCG, invParam.MaxCG, res.n_count);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpMultiSysSolverMR() {}

    Handle< LinearOperator<T> > A;
    MultiSysSolverMRParams invParam;
  };


} // End namespace

#endif 

