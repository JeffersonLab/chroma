// -*- C++ -*-
// $Id: multi_syssolver_mdagm_cg_accumulate.h,v 3.1 2008-09-02 20:10:18 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#ifndef __multi_syssolver_mdagm_cg_accumulate_h__
#define __multi_syssolver_mdagm_cg_accumulate_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/minvcg2_accum.h"


namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverAccumulateCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMMultiSysSolverCGAccumulate : public MdagMMultiSystemSolverAccumulate<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCGAccumulate(Handle< LinearOperator<T> > A_,
			  const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMMultiSysSolverCGAccumulate() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const Real& norm, const multi1d<Real>& residua,
				      const multi1d<Real>& poles, const T& chi) const
      {
	START_CODE();

	/*
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
	  QDPIO::cerr << "MdagMMultiSysSolverCG: shifts incompatible" << endl;
	  QDP_abort(1);
	}

	*/
	SystemSolverResults_t res;
  	MInvCG2Accum(*A, chi, psi, norm, residua, poles, invParam.RsdCG, invParam.MaxCG, res.n_count);


	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMMultiSysSolverCGAccumulate() {}

    Handle< LinearOperator<T> > A;
    SysSolverCGParams invParam;
  };


} // End namespace

#endif 

