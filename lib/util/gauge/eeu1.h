// -*- C++ -*-
// $Id: eeu1.h,v 1.3 2005-01-14 18:42:38 edwards Exp $
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
