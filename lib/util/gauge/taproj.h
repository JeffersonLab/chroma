// -*- C++ -*-
// $Id: taproj.h,v 3.3 2006-08-25 23:56:51 edwards Exp $
/*! \file
 *  \brief Take the traceless antihermitian projection of a color matrix
 */

#ifndef __taproj_h__
#define __taproj_h__

namespace Chroma 
{

  //! Take the traceless antihermitian projection of a color matrix
  /*!
   * \ingroup gauge
   *
   *  a = (1/2)[a - a_dag] - Tr[(1/2)*(a - a_dag)]/Nc
   *
   * that is the anti-hermitian traceless part of a 
   *
   * Arguments:
   *
   *  \param a        LatticeColorMatrix          (Modify)
   */

  void taproj(LatticeColorMatrix& a);

  namespace TaprojEnv {
    extern double getTime();
  }

}  // end namespace Chroma

#endif
