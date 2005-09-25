// $Id: asqtad_cps_wrapper_qprop.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 *  \brief Propagator solver for an even-odd non-preconditioned fermion operator
 *
 *  Solve for the propagator of an even-odd non-preconditioned fermion operator
 */

#ifndef ASQTAD_CPS_WRAPPER_QPROP_H
#define ASQTAD_CPS_WRAPPER_QPROP_H


#include "chromabase.h"
#include "handle.h"
#include "invtype.h"
#include "syssolver.h"
#include "linearop.h"
#include "fermact.h"
#include "state.h"

namespace Chroma
{

  //! QPROP for wrapping CPS asqtad QPROP
  /*! \ingroup qprop
   *
   */
  class AsqtadCPSWrapperQprop : public SystemSolver<LatticeStaggeredFermion>
  {
  public:
    //! Constructor
    /*!
      // Keeping the same interface as for the ordinary staggered 
      // qprop...
      //
      // But the M_ and A_ linop handles are no longer kept
      // (are ignored) -- is there a nice way around this ? 
      // Perhaps not
     */


    AsqtadCPSWrapperQprop(const EvenOddStaggeredTypeFermAct<LatticeStaggeredFermion,multi1d<LatticeColorMatrix> >& S_,
			  Handle< const ConnectState > state, 
			  const InvertParam_t& invParam_);

    //! Destructor is automatic
    ~AsqtadCPSWrapperQprop();

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (LatticeStaggeredFermion& psi, const LatticeStaggeredFermion& chi) const;

  private:
    // Hide default constructor
    AsqtadCPSWrapperQprop() {}

    const Real Mass;
    const InvertParam_t invParam;
    Handle< const ConnectState> state;
  };

}; // End namespace

#endif 

