#ifndef GAUGEACT_IO_H
#define GAUGEACT_IO_H

#include "chromabase.h"
#include "enum_io/enum_io.h"
#include <string>

using namespace QDP;
using namespace std;
using namespace Chroma;

class GaugeActParamsBase {
 public:
  // virtual destructor
  virtual ~GaugeActParamsBase(void) {}

  // Virtual copy function
  virtual GaugeActParamsBase* clone(void) const = 0;

  // get Type
  virtual GaugeActType getType(void) const =0;
};

GaugeActParamsBase* readGaugeActParams(XMLReader& xml, const string& path);
void write(XMLWriter& xml, const string& path, const GaugeActParamsBase& p);



#endif
