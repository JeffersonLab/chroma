// $Id: gaugeacts_aggregate.cc,v 3.7 2006-09-07 04:24:28 edwards Exp $
/*! \file
 *  \brief Generic gauge action wrapper
 */

#include "chromabase.h"

#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"

#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/pg_gaugeact.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/lw_tree_gaugeact.h"
#include "actions/gauge/gaugeacts/lw_1loop_gaugeact.h"
#include "actions/gauge/gaugeacts/rg_gaugeact.h"
#include "actions/gauge/gaugeacts/rbc_gaugeact.h"
#include "actions/gauge/gaugeacts/spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_spectrum_gaugeact.h"

#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"

namespace Chroma
{

  //! Registration aggregator
  namespace GaugeActsEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Register all gauge BCs
      success &= GaugeTypeGaugeBCEnv::registered;

      // Register all gauge states
      success &= CreateGaugeStateEnv::registered;

      // Register gauge actions
      success &= PlaqGaugeActEnv::registered;
      success &= RectGaugeActEnv::registered;
      success &= PgGaugeActEnv::registered;
      success &= WilsonGaugeActEnv::registered;
      success &= LWTreeGaugeActEnv::registered;
      success &= LW1LoopGaugeActEnv::registered;
      success &= RGGaugeActEnv::registered;
      success &= RBCGaugeActEnv::registered;
      success &= SpatialTwoPlaqGaugeActEnv::registered;
      success &= PlaqPlusSpatialTwoPlaqGaugeActEnv::registered;
      success &= AnisoSpectrumGaugeActEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
