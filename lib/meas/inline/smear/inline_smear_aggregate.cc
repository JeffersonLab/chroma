// $Id: inline_smear_aggregate.cc,v 3.2 2006-09-20 20:28:03 edwards Exp $
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
	success &= InlineLinkSmearEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
