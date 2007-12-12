#ifndef UNPREC_S_TPREC_T_WILSTYPE_FERMACT_W_H
#define UNPREC_S_TPREC_T_WILSTYPE_FERMACT_W_H

/*! @file
 * @brief Unpreconditioned spatial Preconditioned Temporal Wilson like fermion action
 */

#include "central_tprec_fermact_w.h"

namespace Chroma {
 //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Spatial, Central Temporal Preconditioned  Wilson-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Spatially Unpreconditioned, Centrally temporally preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class UnprecSpaceCentralPrecTimeWilsonTypeFermAct : public CentralTimePrecFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup
    virtual ~UnprecSpaceCentralPrecTimeWilsonTypeFermAct() {}

    //! Produce a linear operator for this action
    virtual UnprecSpaceCentralPrecTimeLinearOperator<T,P,Q>*  linOp( Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a the gamma5 hermitian version
    virtual LinearOperator<T>*  hermitianLinOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;
  };


}
#endif
