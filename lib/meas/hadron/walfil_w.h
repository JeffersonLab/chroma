// -*- C++ -*-
// $Id: walfil_w.h,v 1.2 2005-01-14 18:42:36 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#ifndef __wallfil_w_h__
#define __wallfil_w_h__

namespace Chroma {

//! Fill a specific color and spin index with 1.0 on a wall
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * \param a            Source fermion
 * \param slice        time slice
 * \param mu           direction of slice
 * \param color_index  Color index
 * \param spin_index   Spin index
 */

void walfil(LatticeFermion& a, int slice, int mu, int color_index, int spin_index);

}  // end namespace Chroma

#endif
