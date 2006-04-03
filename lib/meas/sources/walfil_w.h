// -*- C++ -*-
// $Id: walfil_w.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#ifndef __wallfil_w_h__
#define __wallfil_w_h__

namespace Chroma 
{

  //! Fill a specific color and spin index with 1.0 on a wall
  /*!
   * \ingroup sources
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
