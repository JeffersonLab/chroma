/*! \file
 *  \brief Inline glue measurement aggregator
 */

#include "meas/inline/glue/inline_glue_aggregate.h"
#include "meas/inline/glue/inline_plaquette.h"
#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/glue/inline_qactden.h"
#include "meas/inline/glue/inline_plaq_density.h"
#include "meas/inline/glue/inline_qnaive.h"
#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/glue/inline_fuzwilp.h"
#include "meas/inline/glue/inline_apply_gaugestate.h"
#include "meas/inline/glue/inline_random_transf_gauge.h"
#include "meas/inline/glue/inline_glue_matelem_colorvec.h"
#include "meas/inline/glue/inline_glue_diag_matelem_colorvec.h"
#include "meas/inline/glue/inline_glueball_ops.h"
#include "meas/inline/glue/inline_wilson_flow.h"


#include "meas/inline/glue/inline_gluon_cc.h"
#include "meas/inline/glue/inline_gluon_fraction_OA.h"
#include "meas/inline/glue/inline_gluon_fraction_OB.h"
#include "meas/inline/glue/inline_gluon_npr.h"
#include "meas/inline/glue/inline_gluon_pdf.h"


namespace Chroma
{

  //! Name and registration
  namespace InlineGlueAggregateEnv
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
	success &= InlinePlaquetteEnv::registerAll();
	success &= InlinePolyakovLoopEnv::registerAll();
        success &= InlineQActDenEnv::registerAll();
        success &= InlinePlaqDenEnv::registerAll();
        success &= InlineQTopEnv::registerAll();
	success &= InlineWilsonLoopEnv::registerAll();
	success &= InlineFuzzedWilsonLoopEnv::registerAll();
	success &= InlineRandomTransfGaugeEnv::registerAll();
	success &= InlineGaugeStateEnv::registerAll();
	success &= InlineGaugeStateEnv::registerAll();
	success &= InlineGlueMatElemColorVecEnv::registerAll();
	success &= InlineGlueDiagMatElemColorVecEnv::registerAll();
	success &= InlineGlueballOpsEnv::registerAll();
	success &= InlineWilsonFlowEnv::registerAll();


	success &= InlineGluonMomFracOAEnv::registerAll();
	success &= InlineGluonMomFracOBEnv::registerAll();
	success &= InlineGluonCCEnv::registerAll();
	success &= InlineGluonNPREnv::registerAll();
	success &= InlineGluonPDFEnv::registerAll();
        


	registered = true;
      }
      return success;
    }
  }

}
