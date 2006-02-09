// $Id: seqsource_aggregate_w.cc,v 2.1 2006-02-09 02:25:25 edwards Exp $
/*! \file
 *  \brief All sequential source constructors
 */

#include "meas/hadron/seqsource_aggregate_w.h"

#include "meas/hadron/mesonseqsrc_w.h"
#include "meas/hadron/barseqsrc_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace HadronSeqSourceEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Hadron
      success &= SimpleMesonSeqSourceEnv::registered;
      success &= SimpleBaryonSeqSourceEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
