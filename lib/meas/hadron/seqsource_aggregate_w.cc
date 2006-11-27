// $Id: seqsource_aggregate_w.cc,v 3.2 2006-11-27 04:33:36 edwards Exp $
/*! \file
 *  \brief All sequential source constructors
 */

#include "meas/hadron/seqsource_aggregate_w.h"

#include "meas/hadron/simple_meson_seqsrc_w.h"
#include "meas/hadron/simple_baryon_seqsrc_w.h"
#include "meas/hadron/deriv_meson_seqsrc_w.h"
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
