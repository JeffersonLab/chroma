// -*- C++ -*-
// $Id: syssolver_mdagm_cg_array.h,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#ifndef __syssolver_mdagm_cg_array_h__
#define __syssolver_mdagm_cg_array_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2_array.h"


namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverCGArrayEnv
  {
    //! Name to be used
    extern const std::string name;

    //! Register the syssolver
    extern const bool registered;
  }


  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverCGArray : public MdagMSystemSolverArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverCGArray(Handle< LinearOperatorArray<T> > A_,
			  const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMSysSolverCGArray() {}

    //! Expected length of array index
    int size() const {return A->size();}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const
      {
	START_CODE();

	SystemSolverResults_t res = InvCG2(*A, chi, psi, invParam.RsdCG, invParam.MaxCG);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMSysSolverCGArray() {}

    Handle< LinearOperatorArray<T> > A;
    SysSolverCGParams invParam;
  };

} // End namespace

#endif 

