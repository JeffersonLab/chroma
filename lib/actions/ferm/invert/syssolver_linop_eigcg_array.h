// -*- C++ -*-
// $Id: syssolver_linop_eigcg_array.h,v 1.6 2009-07-17 19:14:46 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system array by EigCG2
 */

#ifndef __syssolver_linop_eigcg_array_h__
#define __syssolver_linop_eigcg_array_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_eigcg_params.h"
#include "actions/ferm/invert/containers.h"

namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGArrayEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2 with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverEigCGArray : public LinOpSystemSolverArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    LinOpSysSolverEigCGArray(Handle< LinearOperatorArray<T> > A_,
			     const SysSolverEigCGParams& invParam_) : 
      MdagM(new MdagMLinOpArray<T>(A_)), A(A_), invParam(invParam_) 
      {
	// NEED to grab the eignvectors from the named buffer here
	if (! TheNamedObjMap::Instance().check(invParam.eigen_id))
	{
	  TheNamedObjMap::Instance().create< LinAlg::RitzPairsArray<T> >(invParam.eigen_id);
	  LinAlg::RitzPairsArray<T>& GoodEvecs = 
	    TheNamedObjMap::Instance().getData< LinAlg::RitzPairsArray<T> >(invParam.eigen_id);

	  if(invParam.Neig_max>0 ){
	    GoodEvecs.init(invParam.Neig_max,MdagM->size());
	  }
	  else{
	    GoodEvecs.init(invParam.Neig,MdagM->size());
	  }
	}
      }

    //! Destructor is automatic
    ~LinOpSysSolverEigCGArray()
      {
	if (invParam.cleanUpEvecs)
	{
	  TheNamedObjMap::Instance().erase(invParam.eigen_id);
	}
      }

    //! Return the subset on which the operator acts
    const Subset& subset(void) const {return A->subset();}

    //! Expected length of array index
    int size(void) const {return A->size(); }

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     *
     * Definitions supplied in the correspond .cc file
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const;

  private:
    // Hide default constructor
    LinOpSysSolverEigCGArray() {}

    Handle< LinearOperatorArray<T> > MdagM;
    Handle< LinearOperatorArray<T> > A;
    SysSolverEigCGParams invParam;
  };

} // End namespace

#endif 

