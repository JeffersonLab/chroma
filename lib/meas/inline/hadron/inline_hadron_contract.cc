// $Id: inline_hadron_contract.cc,v 3.2 2007-06-10 14:49:06 edwards Exp $
/*! \file
 * \brief Inline hadron contraction calculations - for correlators
 *
 * Hadron spectrum calculations. The general version that write output
 * into lime files.
 */

#include "handle.h"
#include "meas/inline/hadron/inline_hadron_contract.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/hadron/hadron_contract.h"
#include "meas/hadron/hadron_contract_factory.h"
#include "meas/hadron/hadron_contract_aggregate.h"
#include "meas/glue/mesplq.h"
#include "util/info/unique_id.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  // Now for all the actual work
  namespace InlineHadronContractEnv 
  { 
    // The callback stuff
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

    const std::string name = "HADRON_CONTRACT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Make sure all hadron contract function objects are registered
	success &= HadronContractEnv::registerAll();

	// Register this object
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //! Reader for parameters
    void read(XMLReader& xml, const string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);
      read(paramtop, "mom_origin", param.mom_origin);
    }


    //! Writer for parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom2_max", param.mom2_max);
      write(xml, "avg_equiv_mom", param.avg_equiv_mom);
      write(xml, "mom_origin", param.mom_origin);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t::Correlator_t& input)
    {
      XMLReader inputtop(xml, path);

      input.xml = readXMLGroup(inputtop, "Correlators", "CorrelatorType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t::Correlator_t& input)
    {
      push(xml, path);

      xml << input.xml.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "output_file", input.output_file);
      read(inputtop, "Correlators", input.correlators);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "output_file", input.output_file);
      write(xml, "Correlators", input.correlators);

      pop(xml);
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

	// Parameters for source construction
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      write(xml_out, "Param", param);
      write(xml_out, "NamedObject", named_obj);
      write(xml_out, "xml_file", xml_file);

      pop(xml_out);
    }





    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "HadronContract");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineHadronContractEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineHadronContractEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "HadronContract");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << " HADRONCONTRACT: Spectroscopy for any type of fermion" << endl;
      QDPIO::cout << "    Volume: " << Layout::lattSize()[0];
      for (int i=1; i<Nd; ++i) {
	QDPIO::cout << " x " << Layout::lattSize()[i];
      }
      QDPIO::cout << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Start the data file output
      XMLBufferWriter file_xml;
      push(file_xml, "HadronContract");
      write(file_xml, "id", uniqueId());  // NOTE: new ID form
      pop(file_xml);

      // Write the scalar data
      QDPFileWriter qio_output(file_xml, params.named_obj.output_file, 
			       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

      // First calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      // Keep an array of all the xml output buffers
      push(xml_out, "HadronMeasurements");

      // Now loop over the various fermion pairs
      for(int lcorr=0; lcorr < params.named_obj.correlators.size(); ++lcorr)
      {
	const GroupXML_t& had_xml = params.named_obj.correlators[lcorr].xml;

	push(xml_out, "elem");

	// Factory construction
	try
	{
	  // Create and use the hadron 2pt object
	  std::istringstream  xml_s(had_xml.xml);
	  XMLReader  hadtop(xml_s);
	
	  Handle<HadronContract> hadronContract(
	    TheHadronContractFactory::Instance().createObject(
	      had_xml.id,
	      hadtop,
	      had_xml.path));

	  // Compute possibly several correlators
	  std::list< Handle<HadronContractResult_t> > hadron_cont =
	    (*hadronContract)(u, "HadronCorrelator", "CorrelatorType");

	  // Run over the output list
	  for(std::list< Handle<HadronContractResult_t> >::const_iterator had_ptr= hadron_cont.begin(); 
	      had_ptr != hadron_cont.end(); 
	      ++had_ptr)
	  {
	    const Handle<HadronContractResult_t>& had_cont = *had_ptr;
	    XMLBufferWriter xml_buf;
	    xml_buf << had_cont->xml;

	    // Save the qio output file entry
	    write(qio_output, xml_buf, had_cont->bin);
	  }
	}
	catch(const std::string& e) 
	{
	  QDPIO::cerr << InlineHadronContractEnv::name << ": Caught Exception in HadronCorrelator: " 
		      << e << endl;
	  QDP_abort(1);
	}

	pop(xml_out);  // array element
      }
      pop(xml_out);  // spectroscopy
      pop(xml_out);  // HadronContract

      // Close data file
      close(qio_output);

      snoop.stop();
      QDPIO::cout << InlineHadronContractEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineHadronContractEnv::name << ": ran successfully" << endl;

      END_CODE();
    }  // func

  } // namespace InlineHadronContractEnv

}  // namespace Chroma
