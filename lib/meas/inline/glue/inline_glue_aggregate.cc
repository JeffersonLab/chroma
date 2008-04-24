/*! \file
 *  \brief Inline glue measurement aggregator
 */

#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/glue/inline_qnaive.h"
#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/glue/inline_fuzwilp.h"
#include "meas/inline/glue/inline_apply_gaugestate.h"
#include "meas/inline/glue/inline_random_transf_gauge.h"

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
        success &= InlineQTopEnv::registerAll();
	success &= InlineWilsonLoopEnv::registerAll();
	success &= InlineFuzzedWilsonLoopEnv::registerAll();
	success &= InlineRandomTransfGaugeEnv::registerAll();
	success &= InlineGaugeStateEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
