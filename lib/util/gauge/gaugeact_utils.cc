#include "chromabase.h"
#include "util/gauge/gaugeact_utils.h"
#include "actions/gauge/wilson_gaugeact.h"

//! Telephone book that creates a GaugeBC from a bas GaugeBCParams class
GaugeAction* getGaugeActFromParams(Handle< GaugeBC > gbc, 
				   const GaugeActParamsBase& b)
{
  switch(b.getType()) {
  case GAUGEACT_WILSON:
    {
      const WilsonGaugeActParams& w_p = 
	dynamic_cast<const WilsonGaugeActParams&>(b);

      return new WilsonGaugeAct(gbc, w_p);
    }
    break;
  default:
    QDPIO::cerr << "Unknown GaugeBCType_t " << b.getType() << endl;
  }
}
