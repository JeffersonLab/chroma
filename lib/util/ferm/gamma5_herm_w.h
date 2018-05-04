// -*- C++ -*-
/*! \file
 *  \brief Gamma5 hermiticity
 */

#ifndef __gamma5_herm_w_h__
#define __gamma5_herm_w_h__

#include "chromabase.h"

namespace Chroma
{
  //! Return gamma_5*adj(source)*gamma_f
  /*!
   * \ingroup ferm
   *
   * \return \f$\gamma_5*source^\dag*\gamma_5\f$
   *
   */
  LatticePropagator gamma5Herm(const LatticePropagator& a);

}


#endif


