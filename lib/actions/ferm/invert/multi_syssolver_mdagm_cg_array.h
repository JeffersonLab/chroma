// -*- C++ -*-
// $Id: multi_syssolver_mdagm_cg_array.h,v 3.4 2008-04-05 19:04:38 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by multi-shift CG
 */

#ifndef __multi_syssolver_mdagm_cg_array_h__
#define __multi_syssolver_mdagm_cg_array_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm.h"
#include "actions/ferm/invert/multi_syssolver_cg_params.h"
#include "actions/ferm/invert/minvcg_array.h"


namespace Chroma
{

  //! CG system solver namespace
  namespace MdagMMultiSysSolverCGArrayEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a CG system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMMultiSysSolverCGArray : public MdagMMultiSystemSolverArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCGArray(Handle< LinearOperatorArray<T> > A_,
			       const MultiSysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMMultiSysSolverCGArray() {}

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
    SystemSolverResults_t operator() (multi1d< multi1d<T> >& psi, 
				      const multi1d<Real>& shifts, 
				      const multi1d<T>& chi) const
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
  	MInvCG(*A, chi, psi, shifts, RsdCG, invParam.MaxCG, res.n_count);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMMultiSysSolverCGArray() {}

    Handle< LinearOperatorArray<T> > A;
    MultiSysSolverCGParams invParam;
  };

} // End namespace

#endif 

