// -*- C++ -*-
// $Id: expm12.h,v 3.1 2006-08-25 23:46:37 edwards Exp $
/*! \file
 *  \brief 12-th order exponentiation of a lattice color matrix
 */

#ifndef __expm12_h__
#define __expm12_h__

namespace Chroma {

//! 12-th order exponentiation of a lattice color matrix
/*!
 * \ingroup gauge
 *
 *  In place  a_ = 1 + a_ + (1/2)*a_^2 + ...+ (1/n!)*(a_)^n  n = power
 *
 * Arguments:
 *
 *  \param a        LatticeColorMatrix          (Modify)
 */

void expm12(LatticeColorMatrix& a);

}  // end namespace Chroma

#endif
