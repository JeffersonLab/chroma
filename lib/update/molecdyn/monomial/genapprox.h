// -*- C++ -*-
/*! @file
 * @brief Wrapper for Remez code
 */

#ifndef __genapprox_h__
#define __genapprox_h__

#include "actions/ferm/fermacts/remez_coeff.h"

namespace Chroma
{
  //! Wrapper for Remez code
  /*! Compute partial fraction expansions for force, action and heat-bath */
  void generateApprox(RemezCoeff_t& fpfe, RemezCoeff_t& spfe, RemezCoeff_t& sipfe,
		      const Real& lower, const Real& upper,
		      int nth_root, int force_degree, int action_degree,
		      int digit_precision);

}
#endif
