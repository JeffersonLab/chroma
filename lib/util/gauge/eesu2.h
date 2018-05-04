// -*- C++ -*-
/*! \file
 *  \brief Exactly exponentiate a SU(2) lie algebra element
 */

#ifndef __eesu2_h__
#define __eesu2_h__

namespace Chroma 
{
  //! Exactly exponentiate a SU(2) lie algebra element
  /*!
   * \ingroup gauge
   *
   *  \param m        LatticeColorMatrix          (Modify)
   */
  void eesu2(LatticeColorMatrix& m);

}  // end namespace Chroma

#endif
