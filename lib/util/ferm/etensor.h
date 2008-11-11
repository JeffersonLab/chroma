// -*- C++ -*-
// $Id: etensor.h,v 3.1 2008-11-11 21:27:42 edwards Exp $
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
   * \return  \f$Q_{\alpha jk} = 0\quad j\ne k, Q_{111}=1/sqrt(2),Q_{122}=-1/sqrt(2), Q_{211}=-1/sqrt(6), Q_{222}=-1/sqrt(6), Q_{233}=2/sqrt(6)\f$
   */
  Real ETensor3d(int alpha, int j, int k);
}  // end namespace Chroma

#endif
