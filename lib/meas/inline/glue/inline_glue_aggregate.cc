// $Id: inline_glue_aggregate.cc,v 1.2 2005-02-10 15:50:47 edwards Exp $
/*! \file
 *  \brief Inline glue measurement aggregator
 */

#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/glue/inline_wilslp.h"

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
      return success;
    }

    const bool registered = registerAll();
  }

}
