#include "chromabase.h"
#include "io/wilson_gaugeact_io.h"

using namespace std;
using namespace QDP;

WilsonGaugeActParams::WilsonGaugeActParams(XMLReader& xml)
{

  GaugeActType_t my_type;
  try {

    read(xml, "./GaugeAct", my_type);
  }
  catch(const string& e ) { 
    QDPIO::cerr << "Caught exception while reading XML : " << e << endl;
    QDP_abort(1);
  }

  if ( my_type != GAUGEACT_WILSON ) { 
    QDPIO::cerr << "Internal error: WilsonGaugeActParams trying to read wrong kind of params: " << my_type << endl;
    QDP_abort(1);
  }

  try { 
    read(xml, "./beta", beta);
  }
  catch(const string& e ) { 
    QDPIO::cerr << "Caught exception while reading XML : " << e << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const WilsonGaugeActParams& p)
{
  try { 
    push(xml, path);
    write(xml, "GaugeAct", p.getType());
    write(xml, "beta", p.getBeta());
    pop(xml);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception while writing XML: " << e << endl;
    QDP_abort(1);
  }
}
