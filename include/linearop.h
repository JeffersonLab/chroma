// -*- C++ -*-
// $Id: linearop.h,v 1.10 2003-04-03 19:28:06 edwards Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

using namespace QDP;

enum LinOpSign {PLUS = 1, MINUS = -1};


//! Linear Operator
/*! @ingroup linop
 *
 * Supports creation and application for linear operators that
 * hold things like Dirac operators, etc.
 */

class LinearOperator
{
public:
  //! Apply the operator onto a source vector
  virtual LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign) const = 0;

  //! Return the subset on which the operator acts
  virtual const Subset& subset() const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~LinearOperator() {}
};


//! Dslash-like Linear Operator
/*! @ingroup linop
 *
 * These are concessions/optimizations for red-black checkboarding 
 */
class DslashLinearOperator : public LinearOperator
{
public:
  //! Apply operator on both checkerboards
  virtual LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign) const
    {
      LatticeFermion d;

      d[rb[0]] = operator()(psi, isign, 0);
      d[rb[1]] = operator()(psi, isign, 1);

      return d;
    }

  virtual LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign, int cb) const = 0;

  //! Virtual destructor to help in cleanup
  virtual ~DslashLinearOperator() {}
};


#endif
