// -*- C++ -*-
// $Id: su3proj.h,v 1.2 2003-12-29 19:49:17 edwards Exp $
/*! \file
 *  \brief Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#ifndef __su3proj_h__
#define __su3proj_h__

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
