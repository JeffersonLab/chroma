// -*- C++ -*-
// $Id: syssolver_linop_OPTeigcg.h,v 1.3 2008-04-01 21:16:18 kostas Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_OPTeigcg_h__
#define __syssolver_linop_OPTeigcg_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_OPTeigcg_params.h"
#include "actions/ferm/invert/containers.h"

namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverOptEigCGEnv
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
  class LinOpSysSolverOptEigCG : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    LinOpSysSolverOptEigCG(Handle< LinearOperator<T> > A_,
			const SysSolverOptEigCGParams& invParam_) : 
      MdagM(new MdagMLinOp<T>(A_)), A(A_), invParam(invParam_) 
      {
	numMatvecs = 0 ;
	// NEED to grab the eignvectors from the named buffer here
	if (! TheNamedObjMap::Instance().check(invParam.eigen_id))
	{
	  TheNamedObjMap::Instance().create< LinAlg::OptEigInfo >(invParam.eigen_id);
	  LinAlg::OptEigInfo& EigInfo = 
	    TheNamedObjMap::Instance().getData< LinAlg::OptEigInfo >(invParam.eigen_id);
	  int N = Layout::sitesOnNode()*Nc*Ns ;
	  int VectorSpaceSize =  Nc*Ns*(A->subset()).numSiteTable();
	  EigInfo.init(invParam.Neig_max, N, VectorSpaceSize) ;
	}
      }

    //! Destructor is automatic
    ~LinOpSysSolverOptEigCG()
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
    LinOpSysSolverOptEigCG() {}
    int numMatvecs ;
    Handle< LinearOperator<T> > MdagM;
    Handle< LinearOperator<T> > A;
    SysSolverOptEigCGParams invParam;
  };

} // End namespace



#endif 

