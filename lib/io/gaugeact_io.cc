#include "chromabase.h"
#include "io/gaugeact_io.h"

#include "io/wilson_gaugeact_io.h"

using namespace QDP;
using namespace std;

GaugeActParamsBase* readGaugeActParams(XMLReader& xml, const string& path)
{
  XMLReader top(xml, path);

  GaugeActType my_type;

  try { 
    read(top, "./GaugeAct", my_type);
  }
  catch(const string& e) {
    QDPIO::cerr << "Caught exception while reading XML : " << e << endl;
    QDP_abort(1);
  }

  switch(my_type) { 
  case GAUGE_ACT_TYPE_WILSON:
    return new WilsonGaugeActParams(top);
    break;
  default:
    QDPIO::cerr << "Unknown GaugeActType " << my_type << endl;
    QDP_abort(1);
  }

  return 0;
}

void write(XMLWriter& xml, const string& path, const GaugeActParamsBase& p) 
{
  GaugeActType my_type = p.getType();

  switch(my_type) {
  case GAUGE_ACT_TYPE_WILSON: 
    {
      const WilsonGaugeActParams& pp = 
	dynamic_cast<const WilsonGaugeActParams&>(p);

      write(xml, path, pp);
    }
    break;
  default:
    QDPIO::cerr << "Unknown GaugeActType " << my_type << endl;
    QDP_abort(1);
  }
}
