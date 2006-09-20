// $Id: seqsource_aggregate_w.cc,v 3.1 2006-09-20 20:28:01 edwards Exp $
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
	// Hadron
	success &= SimpleMesonSeqSourceEnv::registerAll();
	success &= SimpleBaryonSeqSourceEnv::registerAll();
	success &= DerivMesonSeqSourceEnv::registerAll();
	success &= PhotonRhoSeqSourceEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
