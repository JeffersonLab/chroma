// -*- C++ -*-
// $Id: lmdagm_w.h,v 1.5 2003-11-09 22:35:19 edwards Exp $

#ifndef __lmdagm_w_h__
#define __lmdagm_w_h__

#include "linearop.h"

using namespace QDP;

//! M^dag.M linear operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 * Linear operator forming M^dag.M from an operator M
 */

template<typename T>
class lmdagm : public LinearOperator<T>
{
public:
  //! Full constructor
  lmdagm(const LinearOperator<T>& A_) : A(A_) {}

  //! Destructor
  ~lmdagm() {}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A.subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline T operator() (const T& psi, enum PlusMinus isign) const
    {return A(A(psi, PLUS), MINUS);}

private:
  const LinearOperator<T>& A;
};

#endif
