// -*- C++ -*-
// $Id: syssolver_linop_qdp_mg.h, v1.0 2012-04-05 20:45 sdcohen $
/*! \file
 *  \brief Make contact with the QDP clover multigrid solver, transfer
 *         the gauge field, generate the coarse grids, solve systems
 */

#ifndef __syssolver_linop_qdp_mg_h__
#define __syssolver_linop_qdp_mg_h__
#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/qop_mg/syssolver_qop_mg_params.h"


namespace Chroma
{

  //! QDP multigrid system solver namespace
  namespace LinOpSysSolverQOPMGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system using the external QDP multigrid inverter
  /*! \ingroup invert
   */
  template<typename T> // T is the Lattice Fermion type
  class LinOpSysSolverQOPMG : public LinOpSystemSolver<T>
  {
  public:
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;
 
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param state_    The ferm State (Read)
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverQOPMG(Handle< LinearOperator<T> > A_,
			Handle< FermState<T,Q,Q> > state_,
                        const SysSolverQOPMGParams& invParam_);

    //! Destructor finalizes the QDP environment
    ~LinOpSysSolverQOPMG();

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const;


  private:
    // Hide default constructor
    LinOpSysSolverQOPMG() {}
    Handle< FermState<T,Q,Q> > state;
    Handle< LinearOperator<T> > A;
    SysSolverQOPMGParams invParam;
  };

} // End namespace

#endif 

