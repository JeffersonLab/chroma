// $Id: inline_hadron_aggregate.cc,v 3.5 2006-09-20 20:28:02 edwards Exp $
/*! \file
 *  \brief Inline hadron measurement aggregator
 */

#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/hadron/inline_spectrumQll.h"
#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/hadron/inline_make_source_ferm_w.h"
#include "meas/inline/hadron/inline_propagator_w.h"
#include "meas/inline/hadron/inline_propagator_ferm_w.h"
#include "meas/inline/hadron/inline_multi_propagator_w.h"
#include "meas/inline/hadron/inline_seqsource_w.h"
#include "meas/inline/hadron/inline_hadspec_w.h"
#include "meas/inline/hadron/inline_mesonspec_w.h"
#include "meas/inline/hadron/inline_spectrum_w.h"
#include "meas/inline/hadron/inline_spectrumOct_w.h"
#include "meas/inline/hadron/inline_sink_smear_w.h"
#include "meas/inline/hadron/inline_qqq_w.h"
#include "meas/inline/hadron/inline_qqbar_w.h"
#include "meas/inline/hadron/inline_building_blocks_w.h"
#include "meas/inline/hadron/inline_noisy_building_blocks_w.h"
#include "meas/inline/hadron/inline_bar3ptfn_w.h"
//#include "meas/inline/hadron/inline_multipole_w.h"
#include "meas/inline/hadron/inline_mres_w.h"
#include "meas/inline/hadron/inline_qpropqio_w.h"
#include "meas/inline/hadron/inline_qpropadd_w.h"
#include "meas/inline/hadron/inline_qqqNucNuc_w.h"
#include "meas/inline/hadron/inline_stoch_meson_w.h"
#include "meas/inline/hadron/inline_stoch_baryon_w.h"

// Grab all fermacts to make sure they are registered
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineHadronAggregateEnv
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
	// Grab the fermacts
	success &= WilsonTypeFermActsEnv::registerAll();

	// Hadron stuff
	success &= InlineMakeSourceEnv::registerAll();
	success &= InlinePropagatorEnv::registerAll();

	success &= InlineMakeSourceFermEnv::registerAll();
	success &= InlinePropagatorFermEnv::registerAll();

	success &= InlineMultiPropagatorEnv::registerAll();  // save space
	success &= InlineSeqSourceEnv::registerAll();
	success &= InlineHadSpecEnv::registerAll();
	success &= InlineMesonSpecEnv::registerAll();
	success &= InlineSpectrumEnv::registerAll();
	success &= InlineSpectrumOctEnv::registerAll();
	success &= InlineSinkSmearEnv::registerAll();
	success &= InlineQQQEnv::registerAll();
	success &= InlineQQbarEnv::registerAll();
	success &= InlineBuildingBlocksEnv::registerAll();
	success &= InlineNoisyBuildingBlocksEnv::registerAll();
	success &= InlineBar3ptfnEnv::registerAll();
//      success &= InlineMultipoleEnv::registerAll();  // not being used
	success &= InlineMresEnv::registerAll();
	success &= InlineQpropQIOEnv::registerAll();
	success &= InlineQpropAddEnv::registerAll();
	success &= InlineQQQNucNucEnv::registerAll();
	success &= InlineSpectrumQllEnv::registerAll();
	success &= InlineStochMesonEnv::registerAll();
	success &= InlineStochBaryonEnv::registerAll();

	registered = true;
      }
      return success;
    }

  }

}
