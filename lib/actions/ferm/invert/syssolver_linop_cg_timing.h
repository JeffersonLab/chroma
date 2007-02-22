// -*- C++ -*-
// $Id: syssolver_linop_cg_timing.h,v 3.2 2007-02-22 21:11:46 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_cg_timing_h__
#define __syssolver_linop_cg_timing_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2_timing_hacks.h"


namespace Chroma
{

  //! CG system solver namespace
  namespace LinOpSysSolverCGTimingEnv
  {
    //! Name to be used
    extern const std::string name;

    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverCGTiming : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverCGTiming(Handle< LinearOperator<T> > A_,
		 const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverCGTiming() {}

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

	T chi_tmp;
	(*A)(chi_tmp, chi, MINUS);

	SystemSolverResults_t res;  // initialized by a constructor
	for(int i=0; i < invParam.numRestarts; ++i)
	{
	  int n_count = res.n_count;
	  res = InvCG2_timings(*A, chi_tmp, psi, invParam.MaxCG);
	  res.n_count += n_count;
	}

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverCGTiming() {}

    Handle< LinearOperator<T> > A;
    SysSolverCGParams invParam;
  };

} // End namespace

#endif 

