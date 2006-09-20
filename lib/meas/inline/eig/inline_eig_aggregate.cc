// $Id: inline_eig_aggregate.cc,v 3.1 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief Inline eig measurement aggregator
 */

#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/eig/inline_eigbnds.h"
#include "meas/inline/eig/inline_ritz_H_w.h"

// Grab all fermacts to make sure they are registered
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineEigAggregateEnv
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

	// Eig stuff
	success &= InlineEigBndsMdagMEnv::registerAll();
	success &= InlineRitzEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
