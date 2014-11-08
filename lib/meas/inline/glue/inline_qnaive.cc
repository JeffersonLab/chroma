// $Id: inline_qnaive.cc,v 3.4 2009-08-22 20:03:28 edwards Exp $
/*! \file
 *  \brief Inline naive topological charge
 *
 * Author: Christian Hagen
 */

#include "meas/inline/glue/inline_qnaive.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/qnaive.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"


namespace Chroma 
{ 

  namespace InlineQTopEnv 
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

      const std::string name = "QTOP_NAIVE";
    }

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
  

    //! QTop input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	if (paramtop.count("GaugeState") != 0)
	  param.cgs = readXMLGroup(paramtop, "GaugeState", "Name");
	else
	  param.cgs = CreateGaugeStateEnv::nullXMLGroup();
	break;

      default:
	QDPIO::cerr << "Params::Param_t: " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }
    
      read(paramtop, "k5", param.k5);
    }

    //! QTop output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      xml << param.cgs.xml;
      write(xml, "k5", param.k5);

      pop(xml);
    }


    //! QTop input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
    }

    //! QTop output
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

    //NEEDS CHANGES
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      // Grab the object
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "QTop");
      write(xml_out, "update_no", update_no);
      write(xml_out, "k5", params.param.k5);

      Double qtop;
      qtop_naive(u, params.param.k5, qtop);

      write(xml_out, "qtop", qtop);

      pop(xml_out); // pop("QTop");

      END_CODE();
    } 

  }

}
