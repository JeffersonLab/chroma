// -*- C++ -*-
// $Id: walfil_s.h,v 1.1 2003-12-16 02:15:47 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#ifndef __wallfil_s_h__
#define __wallfil_s_h__

//! Fill a specific color and spin index with 1.0 on a wall
/*!
 * \ingroup hadron
 *
 * This routine is specific to Staggered fermions! 
 *
 * Fill a specific color index with 1.0, on sites in a slice
 * where everything has even coordinates.
 *
 * \param a            Source fermion (write)
 * \param slice        time slice
 * \param mu           direction of slice
 * \param color_index  Color index
 */

void walfil(LatticeFermion& a, int slice, int mu, int color_index);

#endif
