// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.1 2007-09-25 21:17:12 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__
#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_eigcg_params.h"
#include "actions/ferm/invert/inv_eigcg2.h"


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
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverEigCG(Handle< LinearOperator<T> > A_,
		 const SysSolverEigCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverEigCG() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
      {
	START_CODE();

	T chi_tmp;
	(*A)(chi_tmp, chi, MINUS);

	SystemSolverResults_t res;  // initialized by a constructor
	{
	  res = InvEigCG2(*A, psi, chi_tmp, eval, evec, 
			  invParam.Neig, invParam.Nmax, 
			  invParam.RsdCG, invParam.MaxCG);

#ifdef CHROMA_DO_ONE_CG_RESTART
	  // Save existing n_count
	  int n_count = res.n_count;

	  // One automatic restart (if enabled)
	  res = InvEigCG2(*A, chi_tmp, psi, invParam.RsdCGRestart, invParam.MaxCGRestart);
	  res.n_count += n_count;
#endif
	
	}

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverEigCG() {}

    Handle< LinearOperator<T> > A;
    SysSolverEigCGParams invParam;
  };

} // End namespace

#endif 

