// -*- C++ -*-
// $Id: linearop.h,v 1.9 2003-04-02 06:59:03 edwards Exp $

/*! @file
 * @brief Linear Operators
 */

#ifndef __linearop_h__
#define __linearop_h__

using namespace QDP;

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
  virtual LatticeFermion operator() (const LatticeFermion& psi, int isign) const;

  //! Return the subset on which the operator acts
  virtual const Subset& subset() const;

  //! Virtual destructor to help with cleanup;
//  virtual ~LinearOperator() {}
};


//! Dslash-like Linear Operator
/*! @ingroup linop
 *
 * These are concessions/optimizations for red-black checkboarding 
 */
class DslashLinearOperator
{
public:
  //! Apply operator on both checkerboards
  virtual LatticeFermion operator() (const LatticeFermion& psi, int isign) const
    {
      LatticeFermion d;

      d[rb[0]] = this->operator()(psi, isign, 0);
      d[rb[1]] = this->operator()(psi, isign, 1);

      return d;
    }

  virtual LatticeFermion operator() (const LatticeFermion& psi, int isign, int cb) const;

  //! Virtual destructor to help in cleanup
//  virtual ~DslashLinearOperator() {}
};


#endif
