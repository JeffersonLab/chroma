/*! \file
 *  \brief Inline polyakov loop
 */

#include "meas/inline/glue/inline_polylp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/polylp.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{ 

  namespace InlinePolyakovLoopEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path)
      {
	Params p(xml_in, path);
	return new InlineMeas(p);
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "POLYAKOV_LOOP";
  
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }

 

  namespace InlinePolyakovLoopEnv 
  { 
    //! PolyakovLoop input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 2:
	if (paramtop.count("GaugeState") != 0)
	  param.cgs = readXMLGroup(paramtop, "GaugeState", "Name");
	else
	  param.cgs = CreateGaugeStateEnv::nullXMLGroup();
	break;

      default:
	QDPIO::cerr << "InlinePolyakovLoopParams::Param_t: " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }
    }

    //! PolyakovLoop output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 2;
      write(xml, "version", version);
      xml << param.cgs.xml;

      pop(xml);
    }


    //! PolyakovLoop input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
    }

    //! PolyakovLoop output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);

      pop(xml);
    }


    // Params
    Params::Params()
    { 
      frequency = 0; 
      param.cgs = CreateGaugeStateEnv::nullXMLGroup();
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Params
	read(paramtop, "Param", param);

	// Ids
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      // Grab the object
      multi1d<LatticeColorMatrix> u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "PolyakovLoop");
      write(xml_out, "update_no", update_no);

      multi1d<DComplex> polyloop;
      polylp(u, polyloop);

      write(xml_out, "poly_loop", polyloop);

      pop(xml_out); // pop("PolyakovLoop");

      END_CODE();
    } 

  }

}
