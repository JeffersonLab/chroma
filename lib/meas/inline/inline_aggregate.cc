// $Id: inline_aggregate.cc,v 1.1 2005-02-07 04:11:27 edwards Exp $
/*! \file
 *  \brief Inline measurement aggregator
 */

#include "meas/inline/inline_aggregate.h"
#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/glue/inline_glue_aggregate.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineGlueAggregateEnv::registered;
      success &= InlineEigAggregateEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
