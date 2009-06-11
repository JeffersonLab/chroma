// -*- C++ -*-
// $Id: syssolver_linop_cg.h,v 3.7 2009-06-11 15:20:54 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_cg_h__
#define __syssolver_linop_cg_h__
#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2.h"


namespace Chroma
{

  //! CG system solver namespace
  namespace LinOpSysSolverCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverCG : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverCG(Handle< LinearOperator<T> > A_,
		 const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverCG() {}

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
	SystemSolverResults_t res;  // initialized by a constructor
	StopWatch swatch;
	swatch.reset();
	swatch.start();

	T chi_tmp;
	(*A)(chi_tmp, chi, MINUS);
	res = InvCG2(*A, chi_tmp, psi, invParam.RsdCG, invParam.MaxCG);

#ifdef CHROMA_DO_ONE_CG_RESTART
	  // Save existing n_count
	int n_count = res.n_count;

	// One automatic restart (if enabled)
	res = InvCG2(*A, chi_tmp, psi, invParam.RsdCGRestart, invParam.MaxCGRestart);
	res.n_count += n_count;
#endif
	swatch.stop();
	double time = swatch.getTimeInSeconds();

	{ 
	  T r;
	  r[A->subset()]=chi;
	  T tmp;
	  (*A)(tmp, psi, PLUS);
	  r[A->subset()] -= tmp;
	  res.resid = sqrt(norm2(r, A->subset()));
	}
	QDPIO::cout << "CG_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
      QDPIO::cout << "CG_SOLVER_TIME: "<<time<< " sec" << endl;

	

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverCG() {}

    Handle< LinearOperator<T> > A;
    SysSolverCGParams invParam;
  };

} // End namespace

#endif 

