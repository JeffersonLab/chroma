// $Id: inline_hadron_aggregate_s.cc,v 3.0 2006-04-03 04:59:03 edwards Exp $
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
    bool registerAll() 
    {
      bool success = true; 

      // Hadron stuff
      success &= InlineStaggeredSpectrumEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
