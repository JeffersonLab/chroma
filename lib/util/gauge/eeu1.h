// -*- C++ -*-
// $Id: eeu1.h,v 3.0 2006-04-03 04:59:12 edwards Exp $
/*! \file
 *  \brief Exactly exponentiate a U(1) lie algebra element
 */

#ifndef __eeu1_h__
#define __eeu1_h__

namespace Chroma 
{

  //! Exactly exponentiate a U(1) lie algebra element
  /*!
   * \ingroup gauge
   *
   *  \param lambda      LatticeColorMatrix          (Modify)
   */
  void eeu1(LatticeColorMatrix& lambda);

}  // end namespace Chroma

#endif
