#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/polylp.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma { 

  namespace InlinePolyakovLoopEnv { 

    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					   const std::string& path) {

      InlinePolyakovLoopParams p(xml_in, path);
      return new InlinePolyakovLoop(p);
    }

    const std::string name = "POLYAKOV_LOOP";

    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

  };

  void 
  InlinePolyakovLoop::operator()(const multi1d<LatticeColorMatrix>& u,
			      const unsigned long update_no,
			      XMLWriter& xml_out) 
  {

    QDPIO::cout << "In POLYLOOP CODE" << endl;
    push(xml_out, "PolyakovLoop");
    write(xml_out, "update_no", update_no);

    multi1d<DComplex> polyloop(Nd);
    for(int mu=0; mu < Nd; mu++) {
      polylp(u, polyloop[mu], mu);
    }

    write(xml_out, "poly_loop", polyloop);

    pop(xml_out); // pop("PolyakovLoop");
    
  } 

};
