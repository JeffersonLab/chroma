// -*- C++ -*-
// $Id: polylinop.h,v 3.0 2006-04-03 04:58:44 edwards Exp $

/*! @file
 * @brief Polynomial Linear Operators
 */

#ifndef __polylinop_h__
#define __polylinop_h__

#include "linearop.h"

namespace Chroma
{

  //----------------------------------------------------------------
  //! Polynomial linear operator including derivatives
  /*! @ingroup linop
   *
   * Support for polynomial linear operators with derivative
   */
  template<typename T, typename P, typename Q>
  class PolyLinearOperator : public DiffLinearOperator<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~PolyLinearOperator() {}

    //! Apply the A or A_dagger piece
    virtual void applyA(T& chi, const T& psi, enum PlusMinus isign) const = 0;
  };


}



#endif
