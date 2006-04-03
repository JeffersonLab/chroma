// $Id: link_smearing_aggregate.cc,v 3.0 2006-04-03 04:59:05 edwards Exp $
/*! \file
 *  \brief All link smearing applicators
 */

#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/ape_link_smearing.h"
#include "meas/smear/hyp_link_smearing.h"
#include "meas/smear/stout_link_smearing.h"

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
      success &= HypLinkSmearingEnv::registered;
      success &= StoutLinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
