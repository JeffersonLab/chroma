// -*- C++ -*-
// $Id: syssolver_linop_eigcg_qdp.h,v 3.1 2008-04-07 04:58:51 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_qdp_h__
#define __syssolver_linop_eigcg_qdp_h__

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
  namespace LinOpSysSolverQDPEigCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2 with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverQDPEigCG : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    LinOpSysSolverQDPEigCG(Handle< LinearOperator<T> > A_,
			   const SysSolverEigCGParams& invParam_) : 
      MdagM(new MdagMLinOp<T>(A_)), A(A_), invParam(invParam_) 
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
    ~LinOpSysSolverQDPEigCG()
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
    LinOpSysSolverQDPEigCG() {}

    Handle< LinearOperator<T> > MdagM;
    Handle< LinearOperator<T> > A;
    SysSolverEigCGParams invParam;
  };

} // End namespace

#endif 

