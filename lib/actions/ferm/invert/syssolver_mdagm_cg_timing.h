// -*- C++ -*-
// $Id: syssolver_mdagm_cg_timing.h,v 3.1 2006-11-16 20:39:16 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#ifndef __syssolver_mdagm_cg_timing_h__
#define __syssolver_mdagm_cg_timing_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2_timing_hacks.h"


namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverCGTimingsEnv
  {
    //! Name to be used
    extern const std::string name;

    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverCGTimings : public MdagMSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverCGTimings(Handle< LinearOperator<T> > A_,
		     const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMSysSolverCGTimings() {}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return A->subset();}

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
	for(int i=0; i < invParam.numRestarts; ++i)
	{
	  int n_count = res.n_count;
	  res = InvCG2_timings(*A, chi, psi, invParam.MaxCG);
	  res.n_count += n_count;
	}

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMSysSolverCGTimings() {}

    Handle< LinearOperator<T> > A;
    SysSolverCGParams invParam;
  };


} // End namespace

#endif 

