// -*- C++ -*-
// $Id: eesu3.h,v 2.0 2005-09-25 21:04:44 edwards Exp $
/*! \file
 *  \brief Exactly exponentiate a SU(3) lie algebra element
 */

#ifndef __eesu3_h__
#define __eesu3_h__

namespace Chroma 
{
  //! Exactly exponentiate a SU(3) lie algebra element
  /*!
   * \ingroup gauge
   *
   *  \param m        LatticeColorMatrix          (Modify)
   */
  void eesu3(LatticeColorMatrix& m);

}  // end namespace Chroma

#endif
