#include "gaugebc_utils.h"

//! Telephone book that creates a GaugeBC from a bas GaugeBCParams class
GaugeBC* getGaugeBCFromParams(const GaugeBCParamsBase& b)
{
  switch(b.getType()) {
  case GAUGEBC_ALL_PERIODIC:
    {
      return new PeriodicGaugeBC;
    }
    break;
  case GAUGEBC_SIMPLE:
    {
      const GaugeBCSimpleParams& p =
	dynamic_cast<const GaugeBCSimpleParams&>(b);
     
      return new SimpleGaugeBC(p);
    }
    break;
  case GAUGEBC_SCHROEDINGER_1LINK:
    { 
      const GaugeBCSchrParams& p = 
	dynamic_cast<const GaugeBCSchrParams&>(b);

      return new Schr1LinkGaugeBC(p);
    }
    break;
  case GAUGEBC_SCHROEDINGER_2LINK:
    {
      const GaugeBCSchrParams& p = 
	dynamic_cast<const GaugeBCSchrParams&>(b);

      return new Schr1LinkGaugeBC(p);
    }
    break;
  default:
    QDPIO::cerr << "Unknown GaugeBCType_t " << b.getType() << endl;
  }

  return 0;
}
