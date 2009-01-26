// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.13 2009-01-26 22:47:05 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by EigCG
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_mdagm.h"

namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by EigCG with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverEigCG : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_          Linear operator ( Read )
     * \param sysSolver_  MdagM system solver ( Read )
     */
    LinOpSysSolverEigCG(Handle< LinearOperator<T> > A_,
			Handle< MdagMSystemSolver<T> > sysSolver_) 
      : A(A_), sysSolver(sysSolver_) {}

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
	T chi_tmp;	
	(*A)(chi_tmp, chi, MINUS);

	return (*sysSolver)(psi, chi_tmp);
      }

  private:

    // Hide default constructor
    Handle< LinearOperator<T> > A;
    Handle< MdagMSystemSolver<T> > sysSolver;
  };
 
} // End namespace



#endif 

