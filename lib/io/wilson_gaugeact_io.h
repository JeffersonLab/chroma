#ifndef WILSON_GAUGEACT_IO_H
#define WILSON_GAUGEACT_IO_H

#include "chromabase.h"
#include "io/gaugeact_io.h"

#include <string>

using namespace QDP;
using namespace std;

class WilsonGaugeActParams : public GaugeActParamsBase {
public: 
  WilsonGaugeActParams(XMLReader& xml);
  ~WilsonGaugeActParams(void) {}
  
  WilsonGaugeActParams(const WilsonGaugeActParams& p) : beta(p.beta) {}

  WilsonGaugeActParams* clone(void) const {
    return new WilsonGaugeActParams(*this);
  }

  GaugeActType getType(void) const { 
    return GAUGE_ACT_TYPE_WILSON;
  }

  const Real& getBeta(void) const { 
    return beta; 
  }
private:
  Real beta;
};

void write(XMLWriter& xml, const string& path, const WilsonGaugeActParams& g);


#endif
