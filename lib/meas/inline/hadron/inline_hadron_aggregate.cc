// $Id: inline_hadron_aggregate.cc,v 3.22 2007-11-16 22:27:33 kostas Exp $
/*! \file
 *  \brief Inline hadron measurement aggregator
 */

#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/hadron/inline_apply_fermstate_w.h"
//#include "meas/inline/hadron/inline_spectrumQll.h"
#include "meas/inline/hadron/inline_static_light_spec_w.h"
#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/hadron/inline_make_source_ferm_w.h"
#include "meas/inline/hadron/inline_propagator_w.h"
#include "meas/inline/hadron/inline_propagator_ferm_w.h"
#include "meas/inline/hadron/inline_multi_propagator_w.h"
#include "meas/inline/hadron/inline_seqsource_w.h"
#include "meas/inline/hadron/inline_seqprop_test_w.h"
#include "meas/inline/hadron/inline_hadspec_w.h"
#include "meas/inline/hadron/inline_mesonspec_w.h"
#include "meas/inline/hadron/inline_hadron_contract.h"
#include "meas/inline/hadron/inline_stag_to_wils.h"
//#include "meas/inline/hadron/inline_spectrum_w.h"
#include "meas/inline/hadron/inline_sink_smear_w.h"
#include "meas/inline/hadron/inline_diquark_w.h"
#include "meas/inline/hadron/inline_qqq_w.h"
#include "meas/inline/hadron/inline_qqq_diquark_w.h"
#include "meas/inline/hadron/inline_qqbar_w.h"
#include "meas/inline/hadron/inline_building_blocks_w.h"
#include "meas/inline/hadron/inline_noisy_building_blocks_w.h"
#include "meas/inline/hadron/inline_bar3ptfn_w.h"
//#include "meas/inline/hadron/inline_multipole_w.h"
#include "meas/inline/hadron/inline_npr_vertex_w.h"
#include "meas/inline/hadron/inline_npr_w.h"
#include "meas/inline/hadron/inline_mres_w.h"
#include "meas/inline/hadron/inline_qpropqio_w.h"
#include "meas/inline/hadron/inline_qpropadd_w.h"
#include "meas/inline/hadron/inline_qqqNucNuc_w.h"
#include "meas/inline/hadron/inline_stoch_meson_w.h"
#include "meas/inline/hadron/inline_stoch_baryon_w.h"
#include "meas/inline/hadron/inline_stoch_group_baryon_w.h"
#include "meas/inline/hadron/inline_gauge_transf_obj.h"

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
	success &= InlineFermStateEnv::registerAll();

	success &= InlineMakeSourceEnv::registerAll();
	success &= InlinePropagatorEnv::registerAll();

	success &= InlineMakeSourceFermEnv::registerAll();
	success &= InlinePropagatorFermEnv::registerAll();

	success &= InlineMultiPropagatorEnv::registerAll();  // save space
	success &= InlineSeqSourceEnv::registerAll();
	success &= InlineSeqPropTestEnv::registerAll();
	success &= InlineHadSpecEnv::registerAll();
	success &= InlineMesonSpecEnv::registerAll();
	success &= InlineHadronContractEnv::registerAll();
//	success &= InlineSpectrumEnv::registerAll();
	success &= InlineStagToWilsEnv::registerAll();
	success &= InlineSinkSmearEnv::registerAll();
	success &= InlineDiquarkEnv::registerAll();
	success &= InlineQQQDiquarkEnv::registerAll();
	success &= InlineQQQEnv::registerAll();
	success &= InlineQQbarEnv::registerAll();
	success &= InlineBuildingBlocksEnv::registerAll();
	success &= InlineNoisyBuildingBlocksEnv::registerAll();
	success &= InlineBar3ptfnEnv::registerAll();
//      success &= InlineMultipoleEnv::registerAll();  // not being used
	success &= InlineNprVertexEnv::registerAll();
	success &= InlineNprEnv::registerAll();
	success &= InlineMresEnv::registerAll();
	success &= InlineGaugeTransfNamedObjEnv::registerAll();
	success &= InlineQpropQIOEnv::registerAll();
	success &= InlineQpropAddEnv::registerAll();
	success &= InlineQQQNucNucEnv::registerAll();
//	success &= InlineSpectrumQllEnv::registerAll();
	success &= InlineStaticLightSpecEnv::registerAll();
	success &= InlineStochMesonEnv::registerAll();
	success &= InlineStochBaryonEnv::registerAll();
	success &= InlineStochGroupBaryonEnv::registerAll();

	registered = true;
      }
      return success;
    }

  }

}
