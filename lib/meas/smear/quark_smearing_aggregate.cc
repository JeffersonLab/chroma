// $Id: quark_smearing_aggregate.cc,v 3.0 2006-04-03 04:59:05 edwards Exp $
/*! \file
 *  \brief All quark smearing
 */

#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/gaus_quark_smearing.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= GausQuarkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
