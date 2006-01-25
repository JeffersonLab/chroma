// $Id: inline_aggregate.cc,v 2.1 2006-01-25 16:50:41 edwards Exp $
/*! \file
 *  \brief Inline measurement aggregator
 */

#include "meas/inline/inline_aggregate.h"
#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/smear/inline_smear_aggregate.h"
#include "meas/inline/io/inline_io_aggregate.h"

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
      success &= InlineSmearAggregateEnv::registered;
      success &= InlineIOAggregateEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
