// -*- C++ -*-
// $Id: multi_syssolver_mdagm_cg_accumulate_array.h,v 3.1 2008-09-06 18:35:35 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by multi-shift CG
 */

#ifndef __multi_syssolver_mdagm_cg_accumulate_array_h__
#define __multi_syssolver_mdagm_cg_accumulate_array_h__

#include "qdp.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate.h"
#include "actions/ferm/invert/multi_syssolver_cg_params.h"
#include "actions/ferm/invert/minvcg_accumulate_array.h"

using namespace QDP;

namespace Chroma
{

  //! CG system solver namespace
  namespace MdagMMultiSysSolverCGAccumulateArrayEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a CG system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMMultiSysSolverCGAccumulateArray : public MdagMMultiSystemSolverAccumulateArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCGAccumulateArray(Handle< LinearOperatorArray<T> > A_,
			       const MultiSysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMMultiSysSolverCGAccumulateArray() {}

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
    SystemSolverResults_t operator() (multi1d<T>& psi, 
				      const Real& norm,
				      const multi1d<Real>& residues,
				      const multi1d<Real>& poles, 
				      const multi1d<T>& chi) const
      {
	START_CODE();


	SystemSolverResults_t res;
	Real RsdCG=invParam.RsdCG[0];

  	MInvCGAccum(*A, chi,psi, norm, residues,poles, RsdCG, invParam.MaxCG, res.n_count);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMMultiSysSolverCGAccumulateArray() {}

    Handle< LinearOperatorArray<T> > A;
    MultiSysSolverCGParams invParam;
  };

} // End namespace

#endif 

