// $Id: inline_link_smear.cc,v 3.5 2007-08-23 19:02:45 edwards Exp $
/*! \file
 *  \brief Inline Link smearing
 */

#include "handle.h"
#include "meas/inline/smear/inline_link_smear.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/gauge/unit_check.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/link_smearing_factory.h"


namespace Chroma 
{ 
  // Reader for out gauge file
  void read(XMLReader& xml, const string& path, InlineLinkSmearEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);
    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "linksmear_id", input.linksmear_id);
  }

  // Reader for out gauge file
  void write(XMLWriter& xml, const string& path, const InlineLinkSmearEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "linksmear_id", input.linksmear_id);

    pop(xml);
  }


  namespace InlineLinkSmearEnv 
  { 
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

    const std::string name = "LINK_SMEAR";

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
      
	// Read program parameters
	link_smearing = readXMLGroup(paramtop, "Param", "LinkSmearingType");

	// Read in the linksmear outfile
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    // Write params
    void
    Params::writeXML(XMLWriter& xml, const std::string& path) 
    {
      push(xml, path);
      
      xml << link_smearing.xml;
      write(xml, "NamedObject", named_obj);

      pop(xml);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Grab the gauge field
      XMLBufferWriter gauge_xml;
      multi1d<LatticeColorMatrix> u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

      push(xml_out, "linksmear");
      write(xml_out, "update_no", update_no);
    
      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables
      MesPlq(xml_out, "Observables", u);

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;
      try 
      {
	std::istringstream  xml_l(params.link_smearing.xml);
	XMLReader  linktop(xml_l);
//	const string link_path = "/Param";
	QDPIO::cout << "Link smearing type = " << params.link_smearing.id << endl;
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,
								       linktop,
								       params.link_smearing.path));
	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << InlineLinkSmearEnv::name << ": error creating link smearing: " << e << endl;
	QDP_abort(1);
      }

      // Check if the smeared gauge field is unitary
      unitarityCheck(u_smr);
  
      // Again calculate some gauge invariant observables
      MesPlq(xml_out, "Link_observables", u_smr);

      // Now store the configuration to a memory object
      {
	XMLBufferWriter file_xml, record_xml;
	push(file_xml, "gauge");
	write(file_xml, "id", int(0));
	pop(file_xml);
	record_xml << gauge_xml;

	// Store the gauge field
	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.linksmear_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.linksmear_id) = u_smr;
	TheNamedObjMap::Instance().get(params.named_obj.linksmear_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.linksmear_id).setRecordXML(record_xml);
      }

      pop(xml_out);


      snoop.stop();
      QDPIO::cout << InlineLinkSmearEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineLinkSmearEnv::name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
