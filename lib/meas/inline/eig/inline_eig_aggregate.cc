// $Id: inline_eig_aggregate.cc,v 1.1 2005-02-07 04:11:27 edwards Exp $
/*! \file
 *  \brief Inline glue measurement aggregator
 */

#include "meas/inline/eig/inline_eig_aggregate.h"
#include "meas/inline/eig/inline_eigbnds.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineEigAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineEigBndsMdagMEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
