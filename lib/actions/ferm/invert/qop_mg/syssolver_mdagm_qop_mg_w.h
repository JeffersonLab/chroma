// -*- C++ -*-
// $Id: syssolver_linop_qdp_mg.h, v1.0 2012-04-05 20:45 sdcohen $
/*! \file
 *  \brief Make contact with the QDP clover multigrid solver, transfer
 *         the gauge field, generate the coarse grids, solve systems
 */

#ifndef __syssolver_mdagm_qdp_mg_h__
#define __syssolver_mdagm_qdp_mg_h__
#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_mdagm.h"
#include "actions/ferm/invert/qop_mg/syssolver_qop_mg_params.h"

#include "actions/ferm/invert/qop_mg/syssolver_linop_qop_mg_w.h"


namespace Chroma
{

  //! QDP multigrid system solver namespace
  namespace MdagMSysSolverQOPMGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system using the external QDP multigrid inverter
  /*! \ingroup invert
   */
  
  class MdagMSysSolverQOPMG : public MdagMSystemSolver<LatticeFermion>
  {
  public:
    typedef LatticeFermion T;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;
 
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param state_    The ferm State (Read)
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverQOPMG(Handle< LinearOperator<T> > A_,
			Handle< FermState<T,Q,Q> > state_,
                        const SysSolverQOPMGParams& invParam_);

    //! Destructor finalizes the QDP environment
    ~MdagMSysSolverQOPMG();

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const;
    SystemSolverResults_t operator()(T& psi, 
					     const T& chi,
					     AbsChronologicalPredictor4D<T>& predictor) const;


  private:
    // Hide default constructor
    MdagMSysSolverQOPMG() {}
    Handle< FermState<T,Q,Q> > state;
    Handle< LinearOperator<T> > A;
    SysSolverQOPMGParams invParam;
    
    //LinOpSysSolverQOPMG<T> Dinv ;
    Handle< LinOpSysSolverQOPMG> Dinv;  // NB: Balint removed <T> template from LinOpSolver as it was causing hassle
    //Handle< LinOpSystemSolver<T> > Dinv;
  };

} // End namespace

#endif 

