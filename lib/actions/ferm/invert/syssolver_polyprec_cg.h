// -*- C++ -*-
// $Id: syssolver_polyprec_cg.h,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Solve a PolyPrec*psi=chi linear system by CG1
 */

#ifndef __syssolver_polyprec_cg_h__
#define __syssolver_polyprec_cg_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_polyprec.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/invcg1.h"


namespace Chroma
{

  //! CG system solver namespace
  namespace PolyPrecSysSolverCGEnv
  {
    //! Name to be used
    extern const std::string name;

    //! Register the syssolver
    extern const bool registered;
  }


  //! Solve a PolyPrec*psi=chi linear system by CG1
  /*! \ingroup invert
   */
  template<typename T>
  class PolyPrecSysSolverCG : public PolyPrecSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    PolyPrecSysSolverCG(Handle< LinearOperator<T> > A_,
			const SysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~PolyPrecSysSolverCG() {}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return A->subset();}

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

	SystemSolverResults_t res = InvCG1(*A, chi_tmp, psi, invParam.RsdCG, invParam.MaxCG);

	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    PolyPrecSysSolverCG() {}

    Handle< LinearOperator<T> > A;
    SysSolverCGParams invParam;
  };

} // End namespace

#endif 

