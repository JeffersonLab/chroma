// -*- C++ -*-
// $Id: expm12.h,v 1.2 2003-12-29 19:52:57 edwards Exp $
/*! \file
 *  \brief 12-th order exponentiation of a lattice color matrix
 */

#ifndef __expm12_h__
#define __expm12_h__

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

#endif
