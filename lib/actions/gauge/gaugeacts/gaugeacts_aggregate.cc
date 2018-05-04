/*! \file
 *  \brief All gauge actions
 */

#include "chromabase.h"

#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"

#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_adjoint_gaugeact.h"
#include "actions/gauge/gaugeacts/character_gaugeact.h"
#include "actions/gauge/gaugeacts/constrained_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_adj_power_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_plaq_power_gaugeact.h"
#include "actions/gauge/gaugeacts/sf_plaq_plus_adjoint_gaugeact.h"
#include "actions/gauge/gaugeacts/sf_plaq_plus_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_plaq_power_gaugeact.h"
#include "actions/gauge/gaugeacts/pg_gaugeact.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/spatial_wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/temporal_wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/wilson_coarse_fine_gaugeact.h"
#include "actions/gauge/gaugeacts/wilson_coarse_fine_1loop_gaugeact.h"
#include "actions/gauge/gaugeacts/lw_tree_gaugeact.h"
#include "actions/gauge/gaugeacts/lw_1loop_gaugeact.h"
#include "actions/gauge/gaugeacts/rg_gaugeact.h"
#include "actions/gauge/gaugeacts/rbc_gaugeact.h"
#include "actions/gauge/gaugeacts/spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_spectrum_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_sym_spatial_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_sym_temporal_gaugeact.h"

#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{

  //! Registration aggregator
  namespace GaugeActsEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Register all gauge BCs
	success &= GaugeTypeGaugeBCEnv::registerAll();

	// Register all gauge states
	success &= CreateGaugeStateEnv::registerAll();

	// Register gauge actions
	success &= PlaqGaugeActEnv::registerAll();
	success &= RectGaugeActEnv::registerAll();
	success &= PgGaugeActEnv::registerAll();
	success &= WilsonGaugeActEnv::registerAll();
	success &= SpatialWilsonGaugeActEnv::registerAll();
	success &= TemporalWilsonGaugeActEnv::registerAll();
	success &= WilsonCoarseFineGaugeActEnv::registerAll();
	success &= WilsonCoarseFine1LoopGaugeActEnv::registerAll();
	success &= LWTreeGaugeActEnv::registerAll();
	success &= LW1LoopGaugeActEnv::registerAll();
	success &= RGGaugeActEnv::registerAll();
	success &= RBCGaugeActEnv::registerAll();
	success &= SpatialTwoPlaqGaugeActEnv::registerAll();
	success &= PlaqPlusSpatialTwoPlaqGaugeActEnv::registerAll();
	success &= PlaqPlusAdjointGaugeActEnv::registerAll();
	success &= CharacterGaugeActEnv::registerAll();
	success &= PlaqPlusTwoPlaqGaugeActEnv::registerAll();
	success &= SFPlaqPlusAdjointGaugeActEnv::registerAll();
	success &= SFPlaqPlusTwoPlaqGaugeActEnv::registerAll();
	success &= AnisoSpectrumGaugeActEnv::registerAll();
	success &= AnisoSymSpatialGaugeActEnv::registerAll();	
	success &= AnisoSymTemporalGaugeActEnv::registerAll();
	success &= PlaqPlusAdjPowerGaugeActEnv::registerAll();
	success &= PlaqPlusPlaqPowerGaugeActEnv::registerAll();
	success &= ConstrainedPlaqGaugeActEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
