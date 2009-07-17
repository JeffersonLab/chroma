// $Id: inline_random_transf_gauge.cc,v 3.2 2009-07-17 21:32:24 edwards Exp $
/*! \file
 *  \brief Do a random gauge transformation on a gauge field
 */

#include "meas/inline/glue/inline_random_transf_gauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "util/gauge/rgauge.h"
#include "meas/glue/mesplq.h"
#include "util/info/proginfo.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  // Reader for out gauge file
  void read(XMLReader& xml, const string& path, InlineRandomTransfGaugeEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "rgauge_id", input.rgauge_id);
    read(inputtop, "gauge_rot_id", input.gauge_rot_id);
  }

  // Reader for out gauge file
  void write(XMLWriter& xml, const string& path, const InlineRandomTransfGaugeEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "rgauge_id", input.rgauge_id);
    write(xml, "gauge_rot_id", input.gauge_rot_id);

    pop(xml);
  }



  namespace InlineRandomTransfGaugeEnv 
  { 
    //! Anonymous namespace
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "RANDOM_GAUGE_TRANSF";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    // Param stuff
    Params::Params()
    { 
      frequency = 0; 
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
      
	// Read in the gfix outfile
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    // Write params
    void
    Params::writeXML(XMLWriter& xml, const std::string& path) 
    {
      push(xml, path);
      
      write(xml, "NamedObject", named_obj);

      pop(xml);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      QDP::StopWatch snoop;
      snoop.reset();
      snoop.start();

      XMLBufferWriter gauge_xml;
      multi1d<LatticeColorMatrix> u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

      push(xml_out, "RandomTransfGauge");
      write(xml_out, "update_no", update_no);
    
      QDPIO::cout << name << ": coulomb gauge fix" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables
      MesPlq(xml_out, "Observables", u);

      // Now do a random gauge transformation fix
      multi1d<LatticeColorMatrix> u_rgauge = u;
      LatticeColorMatrix g;  // the gauge rotation matrices
      rgauge(u_rgauge, g);

      // Calculate observables again. The link is not gauge invariant.
      MesPlq(xml_out, "Random_gauge_transf_observables", u_rgauge);

      // Now store the configuration to a memory object
      {
	XMLBufferWriter file_xml, record_xml;
	push(file_xml, "gauge");
	write(file_xml, "id", int(0));
	pop(file_xml);
	record_xml << gauge_xml;

	// Store the gauge field
	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.rgauge_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.rgauge_id) = u_rgauge;
	TheNamedObjMap::Instance().get(params.named_obj.rgauge_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.rgauge_id).setRecordXML(record_xml);

	// Store the gauge rotation fields
	TheNamedObjMap::Instance().create< LatticeColorMatrix >(params.named_obj.gauge_rot_id);
	TheNamedObjMap::Instance().getData< LatticeColorMatrix >(params.named_obj.gauge_rot_id) = g;
	TheNamedObjMap::Instance().get(params.named_obj.gauge_rot_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_rot_id).setRecordXML(record_xml);
      }

      pop(xml_out);

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
