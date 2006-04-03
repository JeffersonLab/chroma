// -*- C++ -*-
// $Id: etensor.h,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Tensor used for E representations
 */

#ifndef __etensor_h__
#define __etensor_h__

namespace Chroma 
{

  //! Return E antisymmetric tensor
  /*!
   * \ingroup ferm
   *
   * \return  \f$S_{\alpha jk} = 0\quad j\ne k, S_{111}=S_{222}=+1, S_{122}=S_{233}=-1\f$
   */
  int ETensor3d(int alpha, int j, int k);
  
}  // end namespace Chroma

#endif
