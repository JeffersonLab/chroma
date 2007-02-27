#ifndef EO3DPREC_S_TPREC_T_WILSTYPE_FERMACT_W_H
#define EO3DPREC_S_TPREC_T_WILSTYPE_FERMACT_W_H

/*! @file
 * @brief ILUPreconditioned spatial Preconditioned Temporal Wilson like fermion action
 */

#include "wilstype_fermact_w.h"
#include "central_tprec_linop.h" 

namespace Chroma {
  //-------------------------------------------------------------------------------------------
  //! ILUPreconditioned Spatial, Central Temporal Preconditioned  Wilson-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Spatially ILUPreconditioned, Centrally temporally preconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class EO3DPrecSpaceCentralPrecTimeConstDetWilsonTypeFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup
    virtual ~EO3DPrecSpaceCentralPrecTimeConstDetWilsonTypeFermAct() {}

    //! Produce a linear operator for this action
    virtual EO3DPrecSpaceCentralPrecTimeLinearOperator<T,P,Q>*  linOp( Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a the gamma5 hermitian version
    virtual LinearOperator<T>*  hermitianLinOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return quark prop solver, solution of unpreconditioned system
    virtual SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;
  };


}
#endif
