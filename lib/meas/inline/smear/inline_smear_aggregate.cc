// $Id: inline_smear_aggregate.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Inline smear measurement aggregator
 */

#include "meas/inline/smear/inline_smear_aggregate.h"
#include "meas/inline/smear/inline_ape_smear.h"
#include "meas/inline/smear/inline_hyp_smear.h"
#include "meas/inline/smear/inline_stout_smear.h"

#include "meas/inline/smear/inline_link_smear.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineSmearAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineLinkSmearEnv::registered;

      success &= InlineAPESmearEnv::registered;
      success &= InlineHypSmearEnv::registered;
      success &= InlineStoutSmearEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
