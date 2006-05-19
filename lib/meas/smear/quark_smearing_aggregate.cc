// $Id: quark_smearing_aggregate.cc,v 3.1 2006-05-19 15:05:01 edwards Exp $
/*! \file
 *  \brief All quark smearing
 */

#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/no_quark_smearing.h"
#include "meas/smear/gaus_quark_smearing.h"

namespace Chroma
{

  // Registration aggregator
  namespace QuarkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= NoQuarkSmearingEnv::registered;
      success &= GausQuarkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
