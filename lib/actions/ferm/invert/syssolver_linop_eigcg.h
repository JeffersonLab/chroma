// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.6 2007-10-11 19:00:09 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_eigcg_params.h"
#include "actions/ferm/invert/inv_eigcg2.h"

#include "actions/ferm/invert/containers.h"
#include "actions/ferm/invert/lapack_wrapper.h"
#include "actions/ferm/invert/norm_gram_schm.h"


namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGEnv
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
  class LinOpSysSolverEigCG : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    LinOpSysSolverEigCG(Handle< LinearOperator<T> > A_,
			const SysSolverEigCGParams& invParam_) : 
      MdagM(A_), A(A_), invParam(invParam_) 
      {
	// NEED to grab the eignvectors from the named buffer here
	if (! TheNamedObjMap::Instance().check(invParam.eigen_id))
	{
	  TheNamedObjMap::Instance().create< LinAlg::RitzPairs<T> >(invParam.eigen_id);
	  LinAlg::RitzPairs<T>& GoodEvecs = 
	    TheNamedObjMap::Instance().getData< LinAlg::RitzPairs<T> >(invParam.eigen_id);

	  if(invParam.Neig_max>0 ){
	    GoodEvecs.init(invParam.Neig_max);
	  }
	  else{
	    GoodEvecs.init(invParam.Neig);
	  }
	}
      }

    //! Destructor is automatic
    ~LinOpSysSolverEigCG()
      {
	if (invParam.cleanUpEvecs)
	{
	  TheNamedObjMap::Instance().erase(invParam.eigen_id);
	}
      }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     *
     * Definitions supplied in the correspond .cc file
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const;

  private:
    // Hide default constructor
    LinOpSysSolverEigCG() {}

    Handle< LinearOperator<T> > MdagM;
    Handle< LinearOperator<T> > A;
    SysSolverEigCGParams invParam;
  };

} // End namespace

#endif 

