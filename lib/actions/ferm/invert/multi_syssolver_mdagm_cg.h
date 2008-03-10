// -*- C++ -*-
// $Id: multi_syssolver_mdagm_cg.h,v 3.4 2008-03-10 17:32:40 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#ifndef __multi_syssolver_mdagm_cg_h__
#define __multi_syssolver_mdagm_cg_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm.h"
#include "actions/ferm/invert/multi_syssolver_cg_params.h"
#include "actions/ferm/invert/minvcg.h"
#include "actions/ferm/invert/minvcg2.h"


namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGEnv
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
  class MdagMMultiSysSolverCG : public MdagMMultiSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCG(Handle< LinearOperator<T> > A_,
			  const MultiSysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMMultiSysSolverCG() {}

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
	  QDPIO::cerr << "MdagMMultiSysSolverCG: shifts incompatible" << endl;
	  QDP_abort(1);
	}

	SystemSolverResults_t res;
  	MInvCG2(*A, chi, psi, shifts, RsdCG, invParam.MaxCG, res.n_count);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMMultiSysSolverCG() {}

    Handle< LinearOperator<T> > A;
    MultiSysSolverCGParams invParam;
  };


} // End namespace

#endif 

