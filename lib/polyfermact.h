// -*- C++ -*-
// $Id: polyfermact.h,v 1.1 2006-02-10 02:44:55 edwards Exp $

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
  template<typename T, typename P>
  class PolyWilsonTypeFermAct : public WilsonTypeFermAct<T,P>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~PolyWilsonTypeFermAct() {}

    //! Produce a polynomial preconditioned linear operator for this action
    virtual const DiffLinearOperator<T,P>* polyPrecLinOp(Handle<const ConnectState> state) const = 0;

    //! Produce a polynomial linear operator for this action
    virtual const PolyLinearOperator<T,P>* polyLinOp(Handle<const ConnectState> state) const = 0;
  };

}


#endif
