// $Id: inline_aggregate.cc,v 1.2 2005-04-06 04:34:53 edwards Exp $
/*! \file
 *  \brief Inline measurement aggregator
 */

#include "meas/inline/inline_aggregate.h"
#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/hadron/inline_hadron_aggregate.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineEigAggregateEnv::registered;
      success &= InlineGlueAggregateEnv::registered;
      success &= InlineHadronAggregateEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
