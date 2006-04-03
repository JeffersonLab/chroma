// $Id: gaugeacts_aggregate.cc,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Generic gauge action wrapper
 */

#include "chromabase.h"

#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"

#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/pg_gaugeact.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/lw_tree_gaugeact.h"
#include "actions/gauge/gaugeacts/lw_1loop_gaugeact.h"
#include "actions/gauge/gaugeacts/rg_gaugeact.h"
#include "actions/gauge/gaugeacts/rbc_gaugeact.h"


namespace Chroma
{

  //! Registration aggregator
  namespace GaugeActsEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= PlaqGaugeActEnv::registered;
      success &= RectGaugeActEnv::registered;
      success &= PgGaugeActEnv::registered;
      success &= WilsonGaugeActEnv::registered;
      success &= LWTreeGaugeActEnv::registered;
      success &= LW1LoopGaugeActEnv::registered;
      success &= RGGaugeActEnv::registered;
      success &= RBCGaugeActEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
