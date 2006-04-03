// -*- C++ -*-
// $Id: polyfermact.h,v 3.0 2006-04-03 04:58:44 edwards Exp $

/*! @file
 * @brief Class structure for polynomial fermion actions
 */

#ifndef __polyfermact_h__
#define __polyfermact_h__

#include "chromabase.h"
#include "fermact.h"
#include "polylinop.h"

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
  };

}


#endif
