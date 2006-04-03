// -*- C++ -*-
// $Id: z2_src.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Volume source of Z2 noise
 */

#ifndef  Z2_SRC_INC
#define  Z2_SRC_INC 

namespace Chroma 
{
  //! Z2-source
  /*! @ingroup sources */
  void z2_src(LatticeFermion& a);

  //! Z2-source
  /*! @ingroup sources */
  void z2_src(LatticeStaggeredFermion& a);

  //! Z2-source
  /*! @ingroup sources */
  void z2_src(LatticeFermion& a, int slice, int mu);

}  // end namespace Chroma

#endif
