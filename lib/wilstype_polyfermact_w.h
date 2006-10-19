// -*- C++ -*-
// $Id: wilstype_polyfermact_w.h,v 3.1 2006-10-19 16:01:26 edwards Exp $

/*! @file
 * @brief Class structure for polynomial fermion actions
 */

#ifndef __wilstype_polyfermact_w_h__
#define __wilstype_polyfermact_w_h__

#include "wilstype_fermact_w.h"
#include "polylinop.h"
#include "actions/ferm/invert/syssolver_polyprec.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Polynomial Wilson-like fermion actions with derivatives
  /*! @ingroup actions
   *
   * Polynomial like Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class PolyWilsonTypeFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~PolyWilsonTypeFermAct() {}

    //! Produce a polynomial preconditioned linear operator for this action
    virtual DiffLinearOperator<T,P,Q>* polyPrecLinOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a polynomial linear operator for this action
    virtual PolyLinearOperator<T,P,Q>* polyLinOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    virtual PolyPrecSystemSolver<T>* invPolyPrec(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const = 0;
  };

}


#endif
