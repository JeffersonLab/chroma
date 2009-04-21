

#include "meas/inline/pbp/inline_pbp_aggregate.h"
#include "meas/inline/pbp/inline_psibarpsi_w.h"

// Grab all fermacts to make sure they are registered
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{

  namespace InlinePsiBarPsiAggregateEnv
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
	// Grab the fermacts  
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= InlinePsiBarPsiEnv::registerAll();
	
	registered = true;
	  }
	  return success;
	}
	
  }

}