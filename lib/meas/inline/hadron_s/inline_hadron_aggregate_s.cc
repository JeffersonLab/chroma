// $Id: inline_hadron_aggregate_s.cc,v 3.4 2007-05-09 17:19:44 edwards Exp $
/*! \file
 *  \brief Inline hadron measurement aggregator
 */

#include "meas/inline/hadron_s/inline_hadron_aggregate_s.h"
//#include "meas/inline/hadron_s/inline_spectrum_s.h"
#include "meas/inline/hadron_s/inline_make_source_s.h"
#include "meas/inline/hadron_s/inline_propagator_s.h"
#include "meas/inline/hadron_s/inline_sink_smear_s.h"
#include "meas/inline/hadron_s/inline_apply_fermstate_s.h"

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
//	success &= InlineStaggeredSpectrumEnv::registerAll();  // not being used by UKQCD
	success &= InlineStaggeredMakeSourceEnv::registerAll();
	success &= InlineStaggeredPropagatorEnv::registerAll();
	success &= InlineStaggeredSinkSmearEnv::registerAll();
	success &= InlineStaggeredFermStateEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
