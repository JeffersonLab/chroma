// $Id: inline_hadron_aggregate_s.cc,v 2.1 2006-02-26 21:41:20 edwards Exp $
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
