// $Id: inline_hadron_aggregate.cc,v 2.3 2005-10-21 15:55:39 flemingg Exp $
/*! \file
 *  \brief Inline hadron measurement aggregator
 */

#include "meas/inline/hadron/inline_hadron_aggregate.h"
#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/hadron/inline_propagator_w.h"
#include "meas/inline/hadron/inline_multi_propagator_w.h"
#include "meas/inline/hadron/inline_seqsource_w.h"
#include "meas/inline/hadron/inline_hadspec_w.h"
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

// Grab all fermacts to make sure they are registered
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineHadronAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 

      // Grab the fermacts
      success &= WilsonTypeFermActsEnv::registered;

      // Hadron stuff
      success &= InlineMakeSourceEnv::registered;
      success &= InlinePropagatorEnv::registered;
      success &= InlineMultiPropagatorEnv::registered;
      success &= InlineSeqSourceEnv::registered;
      success &= InlineHadSpecEnv::registered;
      success &= InlineSpectrumEnv::registered;
      success &= InlineSpectrumOctEnv::registered;
      success &= InlineSinkSmearEnv::registered;
      success &= InlineQQQEnv::registered;
      success &= InlineQQbarEnv::registered;
      success &= InlineBuildingBlocksEnv::registered;
      success &= InlineNoisyBuildingBlocksEnv::registered;
      success &= InlineBar3ptfnEnv::registered;
//      success &= InlineMultipoleEnv::registered;  // not being used
      success &= InlineMresEnv::registered;
      success &= InlineQpropQIOEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
