#include "chromabase.h"
#include "io/gaugeact_io.h"

#include "io/wilson_gaugeact_io.h"

using namespace QDP;
using namespace std;

void read(XMLReader& xml, const string& path, GaugeActType_t& g)
{
  string token;

  try { 
    read(xml, path, token);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception reading XML: " << e << endl;
    QDP_abort(1);
  }

  if ( token == "WILSON_GAUGE" ) { 
    g = GAUGEACT_WILSON;
  }
  else { 
    QDPIO::cerr << "Unknown gauge action " << token << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const GaugeActType_t& g) 
{
  string token;

  switch(g) { 
  case GAUGEACT_WILSON:
    token = "WILSON_GAUGE" ;
    break;
  default:
    QDPIO::cerr << "Unsupported GaugeActType_t: " << g << endl;
    QDP_abort(1);
  }

  try { 
    write(xml, path, token);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception writing XML " << e << endl;
    QDP_abort(1);
  }
}

GaugeActParamsBase* readGaugeActParams(XMLReader& xml, const string& path)
{
  XMLReader top(xml, path);

  GaugeActType_t my_type;

  try { 
    read(top, "./GaugeAct", my_type);
  }
  catch(const string& e) {
    QDPIO::cerr << "Caught exception while reading XML : " << e << endl;
    QDP_abort(1);
  }

  switch(my_type) { 
  case GAUGEACT_WILSON:
    return new WilsonGaugeActParams(top);
    break;
  default:
    QDPIO::cerr << "Unknown GaugeActType_t " << my_type << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const GaugeActParamsBase& p) 
{
  GaugeActType_t my_type = p.getType();

  switch(my_type) {
  case GAUGEACT_WILSON: 
    {
      const WilsonGaugeActParams& pp = 
	dynamic_cast<const WilsonGaugeActParams&>(p);

      write(xml, path, pp);
    }
    break;
  default:
    QDPIO::cerr << "Unknown GaugeActType_t " << my_type << endl;
    QDP_abort(1);
  }
}
