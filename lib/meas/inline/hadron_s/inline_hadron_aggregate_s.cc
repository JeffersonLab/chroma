// $Id: inline_hadron_aggregate_s.cc,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Inline hadron measurement aggregator
 */

#include "meas/inline/hadron_s/inline_hadron_aggregate_s.h"
#include "meas/inline/hadron_s/inline_spectrum_s.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineStaggeredHadronAggregateEnv
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
	// Hadron stuff
	success &= InlineStaggeredSpectrumEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
