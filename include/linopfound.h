// -*- C++ -*-
// $Id: linopfound.h,v 1.1 2003-04-09 01:20:35 edwards Exp $

/*! @file
 * @brief Foundry class for Linear Operators
 */

#ifndef __linopfound_h__
#define __linopfound_h__

#include "linearop.h"

using namespace QDP;

//! Foundry class for Linear Operators
/*! @ingroup linop
 *
 * Supports creation of various types of linear operators.
 * E.g., a standard operator A, or a version A^dag.A
 */

class LinOpFoundry
{
public:
  //! Create an operator  M
  virtual LinearOperator* linop() const = 0;

  //! Create an operator  M^dag.M
  virtual LinearOperator* lmdagm() const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~LinOpFoundry() {}
};


#endif
