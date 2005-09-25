// -*- C++ -*-
// $Id: remez_coeff.h,v 2.0 2005-09-25 21:04:26 edwards Exp $
/*! \file
 *  \brief Remez algorithm coefficients
 */

#ifndef __remez_coeff_h__
#define __remez_coeff_h__

#include "chromabase.h"

namespace Chroma
{
  //! Convenient structure to package Remez coeffs
  /*! f(x) = norm + \sum_i res[i]/(x+pole[i]) */
  struct RemezCoeff_t
  {
    multi1d<Real>   res;
    multi1d<Real>   pole;
    Real            norm;
  };
}  // namespace Chroma

#endif  // include guard
