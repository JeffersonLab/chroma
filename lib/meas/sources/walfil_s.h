// -*- C++ -*-
// $Id: walfil_s.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#ifndef __wallfil_s_h__
#define __wallfil_s_h__

namespace Chroma {

//! Fill a specific color and spin index with 1.0 on a wall
/*!
 * \ingroup sources
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
 * \param src_index    Index which defines which corner of a cube on
 *                     the source time slice you want your source to
 *                     be on. The mapping from src_index to site is
 *                     lexicographic, i.e: 0 is (0,0,0), 1 is (1,0,0),
 *                     2 is (0,1,0), 3 is (1,1,0), 4 is (0,0,1),
 *                     5 is (1,0,1), 6 is (0,1,1) and 7 is (1,1,1).
 *

 */

void walfil(LatticeStaggeredFermion& a, int slice, int mu, int color_index, int src_index);

}  // end namespace Chroma

#endif
