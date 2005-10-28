// $Id: link_smearing_aggregate.cc,v 1.1 2005-10-28 21:31:04 edwards Exp $
/*! \file
 *  \brief All link smearing constructors
 */

#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/ape_link_smearing.h"

namespace Chroma
{

  //! Registration aggregator
  namespace LinkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= APELinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
