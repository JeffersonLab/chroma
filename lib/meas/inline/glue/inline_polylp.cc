#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/polylp.h"
#include "meas/inline/io/named_objmap.h"


namespace Chroma { 

  namespace InlinePolyakovLoopEnv { 

    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					   const std::string& path)
    {
      InlinePolyakovLoopParams p(xml_in, path);
      return new InlinePolyakovLoop(p);
    }

    const std::string name = "POLYAKOV_LOOP";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };

 
  //! PolyakovLoop input
  void read(XMLReader& xml, const string& path, InlinePolyakovLoopParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! PolyakovLoop output
  void write(XMLWriter& xml, const string& path, const InlinePolyakovLoopParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  // Params
  InlinePolyakovLoopParams::InlinePolyakovLoopParams()
  { 
    frequency = 0; 
  }

  InlinePolyakovLoopParams::InlinePolyakovLoopParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Ids
      read(paramtop, "NamedObject", named_obj);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void 
  InlinePolyakovLoop::operator()(unsigned long update_no,
				 XMLWriter& xml_out) 
  {
    START_CODE();

    // Grab the object
    multi1d<LatticeColorMatrix> u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "PolyakovLoop");
    write(xml_out, "update_no", update_no);

    multi1d<DComplex> polyloop(Nd);
    for(int mu=0; mu < Nd; mu++) {
      polylp(u, polyloop[mu], mu);
    }

    write(xml_out, "poly_loop", polyloop);

    pop(xml_out); // pop("PolyakovLoop");

    END_CODE();
  } 

};
