#ifndef GAUGEACT_IO_H
#define GAUGEACT_IO_H

#include "chromabase.h"
#include <string>

using namespace QDP;
using namespace std;

enum GaugeActType_t { 
  GAUGEACT_WILSON
};

void read(XMLReader& xml, const string&  path, GaugeActType_t& g);
void write(XMLWriter& xml, const string& path, const GaugeActType_t& g);

class GaugeActParamsBase {
 public:
  // virtual destructor
  virtual ~GaugeActParamsBase(void) {}

  // Virtual copy function
  virtual GaugeActParamsBase* clone(void) const = 0;

  // get Type
  virtual const GaugeActType_t getType(void) const =0;
};

GaugeActParamsBase* readGaugeActParams(XMLReader& xml, const string& path);
void write(XMLWriter& xml, const string& path, const GaugeActParamsBase& p);



#endif
