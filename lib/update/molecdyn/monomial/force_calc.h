// -*- C++ -*-
/*! @file
 * @brief Helper function for calculating forces
 */

#ifndef __force_calc_h__
#define __force_calc_h__

#include "chromabase.h"

namespace Chroma
{
  //! Helper function for calculating forces
  /*! @ingroup molecdyn */
  struct ForceCalc_t
  {
    Real   F_sq;     /*!< sum norm2(F) */
    Real   F_avg;    /*!< sum sqrt(norm2(F)) */
    Real   F_max;    /*!< max(localNorm2(F)) */
  };


  //! Helper function for calculating forces
  /*! @ingroup molecdyn */
  ForceCalc_t forceCalc(const multi1d<LatticeColorMatrix>& F);

}
#endif
