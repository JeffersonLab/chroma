// $Id: inline_aggregate.cc,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 *  \brief Inline measurement aggregator
 */

#include "meas/inline/inline_aggregate.h"
#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/gfix/inline_gfix_aggregate.h"
#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/smear/inline_smear_aggregate.h"
#include "meas/inline/io/inline_io_aggregate.h"

#include "meas/inline/hadron_s/inline_hadron_aggregate_s.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineEigAggregateEnv::registered;
      success &= InlineGFixAggregateEnv::registered;
      success &= InlineGlueAggregateEnv::registered;
      success &= InlineHadronAggregateEnv::registered;
      success &= InlineSmearAggregateEnv::registered;
      success &= InlineIOAggregateEnv::registered;

      success &= InlineStaggeredHadronAggregateEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
