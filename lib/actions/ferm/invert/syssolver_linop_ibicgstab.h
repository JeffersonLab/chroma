// -*- C++ -*-
// $Id: syssolver_linop_ibicgstab.h,v 3.2 2009-07-08 18:46:47 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by BICGSTAB
 */

#ifndef __syssolver_linop_ibicgstab_h__
#define __syssolver_linop_ibicgstab_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_bicgstab_params.h"
#include "actions/ferm/invert/invibicgstab.h"

namespace Chroma
{

  //! IBICGSTAB system solver namespace
  namespace LinOpSysSolverIBiCGStabEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by IBICGSTAB
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverIBiCGStab : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverIBiCGStab(Handle< LinearOperator<T> > A_,
		     const SysSolverBiCGStabParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~LinOpSysSolverIBiCGStab() {}

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
      
      StopWatch swatch;

      SystemSolverResults_t res;  // initialized by a constructor
      
      // For now solve with PLUS until we add a way to explicitly
      // ask for MINUS
      swatch.start();

      res = InvIBiCGStab(*A, 
			chi, 
			psi, 
			invParam.RsdBiCGStab, 
			invParam.MaxBiCGStab, 
			PLUS);
      
      
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
      QDPIO::cout << "IBICGSTAB_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << res.resid/sqrt(norm2(chi,A->subset())) << endl;
      QDPIO::cout << "IBICGSTAB_SOLVER_TIME: "<<time<< " sec" << endl;
   
      END_CODE();
      
      return res;
    }


  private:
    // Hide default constructor
    LinOpSysSolverIBiCGStab() {}

    Handle< LinearOperator<T> > A;
    SysSolverBiCGStabParams invParam;
  };

} // End namespace

#endif 

