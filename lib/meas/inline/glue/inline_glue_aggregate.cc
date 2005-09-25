// $Id: inline_glue_aggregate.cc,v 2.0 2005-09-25 21:04:37 edwards Exp $
/*! \file
 *  \brief Inline glue measurement aggregator
 */

#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/glue/inline_fuzwilp.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineGlueAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlinePlaquetteEnv::registered;
      success &= InlinePolyakovLoopEnv::registered;
      success &= InlineWilsonLoopEnv::registered;
      success &= InlineFuzzedWilsonLoopEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
