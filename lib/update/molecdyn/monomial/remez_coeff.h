// -*- C++ -*-
// $Id: remez_coeff.h,v 3.2 2008-05-14 04:13:44 edwards Exp $
/*! \file
 *  \brief Remez algorithm coefficients
 */

#ifndef __remez_coeff_h__
#define __remez_coeff_h__

#include "chromabase.h"

namespace Chroma
{
  //! Convenient structure to package Remez coeffs
  /*! @ingroup monomial
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
