// -*- C++ -*-
// $Id: syssolver_linop_cg_array.h,v 3.6 2008-04-05 19:04:38 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_cg_array_h__
#define __syssolver_linop_cg_array_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2_array.h"


namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverCGArrayEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverCGArray : public LinOpSystemSolverArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverCGArray(Handle< LinearOperatorArray<T> > A_,
			  const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverCGArray() {}

    //! Expected length of array index
    int size() const {return A->size();}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const
      {
	START_CODE();

	multi1d<T> chi_tmp(size());
	(*A)(chi_tmp, chi, MINUS);

	SystemSolverResults_t res;  // initialized by a constructor
	{

	  res = InvCG2(*A, chi_tmp, psi, invParam.RsdCG, invParam.MaxCG);
	

#ifdef CHROMA_DO_ONE_CG_RESTART

	  int n_count = res.n_count;

	  // One automatic restart (if enabled)
	  res = InvCG2(*A, chi_tmp, psi, invParam.RsdCGRestart, invParam.MaxCGRestart);
	  res.n_count += n_count;
#endif

	}

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverCGArray() {}

    Handle< LinearOperatorArray<T> > A;
    SysSolverCGParams invParam;
  };


} // End namespace

#endif 

