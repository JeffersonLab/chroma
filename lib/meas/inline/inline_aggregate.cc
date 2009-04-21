// $Id: inline_aggregate.cc,v 3.3 2009-04-21 01:21:41 eneil Exp $
/*! \file
 *  \brief Inline measurement aggregator
 */

#include "meas/inline/inline_aggregate.h"
#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/gfix/inline_gfix_aggregate.h"
#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/schrfun/inline_schrfun_aggregate.h"
#include "meas/inline/smear/inline_smear_aggregate.h"
#include "meas/inline/io/inline_io_aggregate.h"
#include "meas/inline/pbp/inline_pbp_aggregate.h"

#include "meas/inline/hadron_s/inline_hadron_aggregate_s.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineAggregateEnv
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
	success &= InlineEigAggregateEnv::registerAll();
	success &= InlineGFixAggregateEnv::registerAll();
	success &= InlineGlueAggregateEnv::registerAll();
	success &= InlineHadronAggregateEnv::registerAll();
	success &= InlineSchrFunAggregateEnv::registerAll();
	success &= InlineSmearAggregateEnv::registerAll();
	success &= InlineIOAggregateEnv::registerAll();
	success &= InlinePsiBarPsiAggregateEnv::registerAll();

	success &= InlineStaggeredHadronAggregateEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
