// $Id: inline_smear_aggregate.cc,v 3.1 2006-08-11 18:11:59 edwards Exp $
/*! \file
 *  \brief Inline smear measurement aggregator
 */

#include "meas/inline/smear/inline_smear_aggregate.h"
#include "meas/inline/smear/inline_link_smear.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineSmearAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineLinkSmearEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
