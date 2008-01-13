// -*- C++ -*-
// $Id: syssolver_linop_eigcg_array.h,v 1.2 2008-01-13 22:43:54 edwards Exp $
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
    //! Name to be used
    extern const std::string name;

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
	  TheNamedObjMap::Instance().create< LinAlg::RitzPairs< multi1d<T> > >(invParam.eigen_id);
	  LinAlg::RitzPairs< multi1d<T> >& GoodEvecs = 
	    TheNamedObjMap::Instance().getData< LinAlg::RitzPairs< multi1d<T> > >(invParam.eigen_id);

	  if(invParam.Neig_max>0 ){
	    GoodEvecs.init(invParam.Neig_max);
	  }
	  else{
	    GoodEvecs.init(invParam.Neig);
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
    const Subset& subset() const {return A->subset();}

    //! Expected length of array index
    int size() const {return A->size();}

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

