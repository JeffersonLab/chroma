// -*- C++ -*-
// $Id: su3proj.h,v 1.1 2003-03-28 03:53:39 edwards Exp $

#ifndef __su3proj__
#define __su3proj__

//! Project a GL(3,C) color matrix onto SU(3)
/*!
 * Arguments:
 *
 *  \param u            the projected SU(3) Matrix (Modify)
 *  \param w            matrix against which to maximize (Read)
 *  \param su2_index    SU(2) subgroup index (Read)
 */

void su3proj(LatticeColorMatrix& u, const LatticeColorMatrix& w, int su2_index);

#endif
