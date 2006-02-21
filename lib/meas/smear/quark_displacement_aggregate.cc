// $Id: quark_displacement_aggregate.cc,v 1.1 2006-02-21 06:45:40 edwards Exp $
/*! \file
 *  \brief All quark displacement
 */

#include "meas/smear/quark_displacement_aggregate.h"

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

      success &= SimpleQuarkDisplacementEnv::registered;
      success &= DerivQuarkDisplacementEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
