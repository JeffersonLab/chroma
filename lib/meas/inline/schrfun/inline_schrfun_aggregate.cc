// $Id: inline_schrfun_aggregate.cc,v 1.2 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Inline Schroedinger functional measurement aggregator
 */

#include "meas/inline/schrfun/inline_schrfun_aggregate.h"

#include "meas/inline/schrfun/inline_sfpcac_w.h"

// Grab all fermacts to make sure they are registered
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineSchrFunAggregateEnv
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
	
	// Schrfun stuff
	success &= InlineSFpcacEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
