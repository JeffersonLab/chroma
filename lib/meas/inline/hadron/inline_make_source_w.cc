// $Id: inline_make_source_w.cc,v 3.6 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of make_source
 *
 * Construct source for propagator calculations
 */

#include "handle.h"
#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/sources/source_const_factory.h"
#include "meas/sources/source_const_aggregate.h"

#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineMakeSourceEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMakeSource(InlineMakeSourceParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MAKE_SOURCE";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= QuarkSourceConstructionEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }


  //! MakeSource input
  void read(XMLReader& xml, const string& path, InlineMakeSourceParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
  }

  //! MakeSource output
  void write(XMLWriter& xml, const string& path, const InlineMakeSourceParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "source_id", input.source_id);

    pop(xml);
  }


  // Param stuff
  InlineMakeSourceParams::InlineMakeSourceParams() { frequency = 0; }

  InlineMakeSourceParams::InlineMakeSourceParams(XMLReader& xml_in, const std::string& path) 
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

      // Named object output location
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


  // Writer
  void
  InlineMakeSourceParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the buffer ids
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  // Function call
  void 
  InlineMakeSource::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "make_source");
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
  InlineMakeSource::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlineMakeSourceEnv::name << ": propagator source constructor" << endl;

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    multi1d<LatticeColorMatrix> u;
    try
    {
      u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": error extracting gauge field: " << e 
		  << endl;
      QDP_abort(1);
    }

    // Save the initial state of the RNG
    QDP::Seed ran_seed;
    QDP::RNG::savern(ran_seed);

    push(xml_out, "make_source");
    write(xml_out, "update_no", update_no);

    proginfo(xml_out);    // Print out basic program info

    // Current state of the seed
    write(xml_out, "RNG", ran_seed);

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //
    // Initialize source
    //
    LatticePropagator quark_source;

    try
    {
      std::istringstream  xml_s(params.param.source.xml);
      XMLReader  sourcetop(xml_s);
      QDPIO::cout << "Source = " << params.param.source.id << endl;

      Handle< QuarkSourceConstruction<LatticePropagator> >
	sourceConstruction(ThePropSourceConstructionFactory::Instance().createObject(params.param.source.id,
										     sourcetop,
										     params.param.source.path));
      quark_source = (*sourceConstruction)(u);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": Caught Exception creating source: " << e << endl;
      QDP_abort(1);
    }


    // Sanity check - write out the norm2 of the source in the Nd-1 direction.
    // Use this for any possible verification.
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_source),
					     phases.getSet());

      push(xml_out, "Source_correlator");
      write(xml_out, "source_corr", source_corr);
      pop(xml_out);
    }
 

    // Now write the source
    try
    {
      QDPIO::cout << "Attempt to update source" << endl;

      XMLBufferWriter file_xml;
      push(file_xml, "make_source");
      write(file_xml, "id", uniqueId());  // NOTE: new ID form
      pop(file_xml);

      XMLBufferWriter record_xml;

      // Construct the appropriate header
      {
	MakeSourceProp_t prop_header;
	prop_header.source_header = params.param;
	prop_header.gauge_header  = gauge_xml.printCurrentContext();
	write(record_xml, "MakeSource", prop_header);
      }
    
      // Store the source
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.source_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id) = quark_source;
      TheNamedObjMap::Instance().get(params.named_obj.source_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).setRecordXML(record_xml);

      QDPIO::cout << "Source successfully update" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": error message: " << e << endl;
      QDP_abort(1);
    }
    
    pop(xml_out);  // make_source

//    // Reset the seed
//    QDP::RNG::setrn(ran_seed);

    snoop.stop();
    QDPIO::cout << InlineMakeSourceEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineMakeSourceEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
