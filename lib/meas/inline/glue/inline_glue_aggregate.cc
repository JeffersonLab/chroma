// $Id: inline_glue_aggregate.cc,v 1.1 2005-02-07 04:11:27 edwards Exp $
/*! \file
 *  \brief Inline glue measurement aggregator
 */

#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/glue/inline_polylp.h"

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
      return success;
    }

    const bool registered = registerAll();
  }

}
