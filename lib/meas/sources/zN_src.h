// -*- C++ -*-
/*! \file
 *  \brief Volume source of Z(N) noise
 */

#ifndef  ZN_SRC_INC
#define  ZN_SRC_INC 

namespace Chroma 
{
  //! Z(N)-rng
  /*! @ingroup sources */
  Complex zN_rng(int N);

  //! Z(N)-source
  /*! @ingroup sources */
  void zN_src(LatticeFermion& a, int N);

}  // end namespace Chroma

#endif
