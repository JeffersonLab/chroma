// -*- C++ -*-
/*! \file
 *  \brief A fake inverter
 */

#ifndef __syssolver_linop_fake_h__
#define __syssolver_linop_fake_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"


namespace Chroma
{

  //! CG system solver namespace
  namespace LinOpSysSolverFakeEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverFake : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     */
    LinOpSysSolverFake(Handle< LinearOperator<T> > A_) : A(A_)
      {}

    //! Destructor is automatic
    ~LinOpSysSolverFake() {}

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

        psi[A->subset()] = chi;

        res.n_count = 0;
        res.resid   = 0;

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    LinOpSysSolverFake() {}

    Handle< LinearOperator<T> > A;
  };

} // End namespace

#endif 

