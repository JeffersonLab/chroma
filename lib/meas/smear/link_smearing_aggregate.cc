// $Id: link_smearing_aggregate.cc,v 1.2 2005-11-07 06:40:55 edwards Exp $
/*! \file
 *  \brief All link smearing applicators
 */

#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/ape_link_smearing.h"

namespace Chroma
{

  // Registration aggregator
  namespace LinkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // link smearing
      success &= APELinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
