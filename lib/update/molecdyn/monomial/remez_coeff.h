// -*- C++ -*-
// $Id: remez_coeff.h,v 3.1 2007-04-17 03:13:04 edwards Exp $
/*! \file
 *  \brief Remez algorithm coefficients
 */

#ifndef __remez_coeff_h__
#define __remez_coeff_h__

#include "chromabase.h"

namespace Chroma
{
  //! Convenient structure to package Remez coeffs
  /*! @ingroup molecdyn
   *
   * f(x) = norm + \sum_i res[i]/(x+pole[i]) 
   */
  struct RemezCoeff_t
  {
    multi1d<Real>   res;
    multi1d<Real>   pole;
    Real            norm;
  };
}  // namespace Chroma

#endif  // include guard
