// -*- C++ -*-
/*! \file
 *  \brief Basis rotation matrix from Dirac to Degrand-Rossi (and reverse)
 */

#ifndef __tdiractodr_h__
#define __tdiractodr_h__

#include "chromabase.h"

namespace Chroma 
{

  //! The Dirac to Degrand-Rossi spin transformation matrix
  /*!
   * \ingroup ferm
   *
   * Return the similarity transformation matrix from 
   * Euclidean Dirac to Euclidean Degrand-Rossi basis
   *
   * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
   */

  SpinMatrixD TDiracToDRMat();

}  // end namespace Chroma

#endif
