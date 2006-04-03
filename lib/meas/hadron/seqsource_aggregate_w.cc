// $Id: seqsource_aggregate_w.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief All sequential source constructors
 */

#include "meas/hadron/seqsource_aggregate_w.h"

#include "meas/hadron/mesonseqsrc_w.h"
#include "meas/hadron/barseqsrc_w.h"
#include "meas/hadron/derivmesonseqsrc_w.h"
#include "meas/hadron/photon_seqsrc_w.h"

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
      success &= DerivMesonSeqSourceEnv::registered;
      success &= PhotonRhoSeqSourceEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
