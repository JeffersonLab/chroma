// -*- C++ -*-
// $Id: unprec_wilstype_fermact_w.h,v 3.1 2006-10-19 16:01:26 edwards Exp $

/*! @file
 * @brief Wilson-like fermion actions
 */

#ifndef __unprec_wilstype_fermact_w_h__
#define __unprec_wilstype_fermact_w_h__

#include "wilstype_fermact_w.h"

namespace Chroma
{

  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class UnprecWilsonTypeFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermAct() {}

    //! Produce a linear operator for this action
    virtual UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Unpreconditioned Wilson-like fermion actions in extra dims with derivatives
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   * Here, use arrays of matter fields.
   */
  template<typename T, typename P, typename Q>
  class UnprecWilsonTypeFermAct5D : public WilsonTypeFermAct5D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~UnprecWilsonTypeFermAct5D() {}

    //! Produce a linear operator for this action
    virtual UnprecLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a Pauli-Villars linear operator for this action
    virtual UnprecLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const = 0;
  };

}


#endif
