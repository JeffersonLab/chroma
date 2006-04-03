// -*- C++ -*-
// $Id: antisymtensor.h,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Compute anti-symmetric tensors
 */

#ifndef __antisymtensor_h__
#define __antisymtensor_h__

namespace Chroma 
{

  //! Return 3d antisymmetric tensor
  /*!
   * \ingroup ferm
   *
   * \return  \f$\epsilon_{ijk}\f$
   */
  int antiSymTensor3d(int i, int j, int k);
  
}  // end namespace Chroma

#endif
