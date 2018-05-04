/*! \file
 *  \brief Inline gauge fixing measurement aggregator
 */

#include "meas/inline/gfix/inline_gfix_aggregate.h"
#include "meas/inline/gfix/inline_coulgauge.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineGFixAggregateEnv
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
	success &= InlineCoulGaugeEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
