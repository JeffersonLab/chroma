// $Id: quark_displacement_aggregate.cc,v 3.0 2006-04-03 04:59:05 edwards Exp $
/*! \file
 *  \brief All quark displacements
 */

#include "meas/smear/quark_displacement_aggregate.h"

#include "meas/smear/no_quark_displacement.h"
#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/deriv_quark_displacement_w.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkDisplacementEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= NoQuarkDisplacementEnv::registered;
      success &= SimpleQuarkDisplacementEnv::registered;
      success &= DerivQuarkDisplacementEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
