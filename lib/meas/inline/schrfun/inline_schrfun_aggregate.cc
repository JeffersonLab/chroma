// $Id: inline_schrfun_aggregate.cc,v 1.1 2006-04-10 21:17:05 edwards Exp $
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
    bool registerAll() 
    {
      bool success = true; 

      // Grab the fermacts
      success &= WilsonTypeFermActsEnv::registered;

      // Schrfun stuff
      success &= InlineSFpcacEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
