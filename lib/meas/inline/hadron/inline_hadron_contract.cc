// $Id: inline_hadron_contract.cc,v 3.5 2007-06-12 17:55:44 edwards Exp $
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


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "output_file", input.output_file);
      QDPIO::cout << "Read contraction list" << endl;
      input.correlators = readXMLArrayGroup(inputtop, "Contractions", "ContractionType");
      QDPIO::cout << "Finished reading correlators" << endl;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "output_file", input.output_file);

      for(int i=0; i < input.correlators.size(); ++i)
      {
	xml << input.correlators[i].xml;
      }

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

      QDPIO::cout << InlineHadronContractEnv::name 
		  << ": hadron contracts and spectroscopy for any type of fermion" 
		  << endl;
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
	const GroupXML_t& had_xml = params.named_obj.correlators[lcorr];

	push(xml_out, "elem");
	write(xml_out, "Input", had_xml.xml);

	// Factory construction
	try
	{
//	  QDPIO::cout << "xml input = XX" << had_xml.xml << "XX" << endl;
	  QDPIO::cout << "Contractions for id = " << had_xml.id << endl;

	  // Create and use the hadron 2pt object
	  std::istringstream  xml_s(had_xml.xml);
	  XMLReader  hadtop(xml_s);
	
	  Handle<HadronContract> hadronContract(
	    TheHadronContractFactory::Instance().createObject(
	      had_xml.id,
	      hadtop,
	      had_xml.path));

	  // Compute possibly several correlators
	  QDPIO::cout << InlineHadronContractEnv::name << ": start list" << endl;

	  std::list< Handle<HadronContractResult_t> > hadron_cont =
	    (*hadronContract)(u, "HadronContraction", "ContractionType");

	  QDPIO::cout << InlineHadronContractEnv::name << ": finished list" << endl;

	  push(xml_out, "HadronContractions");

	  // Run over the output list
	  for(std::list< Handle<HadronContractResult_t> >::const_iterator had_ptr= hadron_cont.begin(); 
	      had_ptr != hadron_cont.end(); 
	      ++had_ptr)
	  {
	    const Handle<HadronContractResult_t>& had_cont = *had_ptr;

	    push(xml_out, "elem");

	    // Save regression output in output file 
	    // NOTE: what is under the "Diagnostic" tag is solely up
	    // to the user. You can jam whatever you want into
	    // here. It is used for the regressions to latch onto
	    // something from the output since otherwise it is all in binary.
	    // The input is written as well to give some context
	    write(xml_out, "RecordXML", had_cont->xml);
	    write(xml_out, "Diagnostic", had_cont->xml_regres);
	    
	    // Save the qio output file entry
	    write(qio_output, had_cont->xml, had_cont->bin);

	    pop(xml_out);  // array element
	  }

	  pop(xml_out);  // HadronContractions
	}
	catch(const std::string& e) 
	{
	  QDPIO::cerr << InlineHadronContractEnv::name << ": Caught Exception in HadronCorrelator: " 
		      << e << endl;
	  QDP_abort(1);
	}

	pop(xml_out);  // array element
      }
      pop(xml_out);  // HadronMeasurements
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
