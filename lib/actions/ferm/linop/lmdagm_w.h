// -*- C++ -*-
// $Id: lmdagm_w.h,v 1.2 2003-08-09 04:18:47 edwards Exp $

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

class lmdagm : public LinearOperator
{
public:
  //! Full constructor
  lmdagm(const LinearOperator& _A) : A(_A) {}

  //! Destructor
  ~lmdagm() {}

  //! Subset comes from underlying operator
  inline const OrderedSubset& subset() const {return A.subset();}

  //! Apply the operator onto a source vector
  /*! For this operator, the sign is ignored */
  inline LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign) const
    {return A(A(psi, PLUS), MINUS);}

private:
  const LinearOperator& A;
};

#endif
