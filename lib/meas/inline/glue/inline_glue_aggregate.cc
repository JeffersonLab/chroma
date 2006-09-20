// $Id: inline_glue_aggregate.cc,v 3.1 2006-09-20 20:28:01 edwards Exp $
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
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= InlinePlaquetteEnv::registerAll();
	success &= InlinePolyakovLoopEnv::registerAll();
	success &= InlineWilsonLoopEnv::registerAll();
	success &= InlineFuzzedWilsonLoopEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
