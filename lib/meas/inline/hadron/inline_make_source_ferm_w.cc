// $Id: inline_make_source_ferm_w.cc,v 3.4 2007-02-25 22:39:28 edwards Exp $
/*! \file
 * \brief Inline construction of make_source for lattice fermions
 *
 * Construct source for propagator calculations
 */

#include "handle.h"
#include "meas/inline/hadron/inline_make_source_ferm_w.h"
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
  namespace InlineMakeSourceFermEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMakeSourceFerm(InlineMakeSourceFermParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MAKE_SOURCE_FERM";

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


  //! MakeSourceFerm input
  void read(XMLReader& xml, const string& path, InlineMakeSourceFermParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
  }

  //! MakeSourceFerm output
  void write(XMLWriter& xml, const string& path, const InlineMakeSourceFermParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "source_id", input.source_id);

    pop(xml);
  }


  // Param stuff
  InlineMakeSourceFermParams::InlineMakeSourceFermParams() { frequency = 0; }

  InlineMakeSourceFermParams::InlineMakeSourceFermParams(XMLReader& xml_in, const std::string& path) 
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


  //! Make source input
  void read(XMLReader& xml, const string& path, InlineMakeSourceFermParams& input)
  {
    InlineMakeSourceFermParams tmp(xml, path);
    input = tmp;
  }

  //! Make source output
  void write(XMLWriter& xml, const string& path, const InlineMakeSourceFermParams& input)
  {
    push(xml, path);
    
    // Parameters for source construction
    Chroma::write(xml, "Param", input.param);

    // Write out the buffer ids
    Chroma::write(xml, "NamedObject", input.named_obj);

    pop(xml);
  }


  // Function call
  void 
  InlineMakeSourceFerm::operator()(unsigned long update_no,
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
  InlineMakeSourceFerm::func(unsigned long update_no,
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
      QDPIO::cerr << InlineMakeSourceFermEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMakeSourceFermEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    // Save the initial state of the RNG
    QDP::Seed ran_seed;
    QDP::RNG::savern(ran_seed);

    push(xml_out, "make_source");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineMakeSourceFermEnv::name << ": propagator source constructor" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Current state of the seed
    write(xml_out, "RNG", ran_seed);

    // Write out the input
    write(xml_out, "Input", params);

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //
    // Initialize source
    //
    LatticeFermion quark_source;

    try
    {
      std::istringstream  xml_s(params.param.source.xml);
      XMLReader  sourcetop(xml_s);
      QDPIO::cout << "Source = " << params.param.source.id << endl;

      Handle< QuarkSourceConstruction<LatticeFermion> >
	sourceConstruction(TheFermSourceConstructionFactory::Instance().createObject(params.param.source.id,
										     sourcetop,
										     params.param.source.path));
      quark_source = (*sourceConstruction)(u);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << InlineMakeSourceFermEnv::name << ": Caught Exception creating source: " << e << endl;
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
      push(record_xml, "MakeSource");
      write(record_xml, "PropSource", params.param);
      write(record_xml, "RNG", ran_seed);
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);
    
      // Store the source
      TheNamedObjMap::Instance().create<LatticeFermion>(params.named_obj.source_id);
      TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.source_id) = quark_source;
      TheNamedObjMap::Instance().get(params.named_obj.source_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).setRecordXML(record_xml);

      QDPIO::cout << "Source successfully update" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineMakeSourceFermEnv::name << ": dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMakeSourceFermEnv::name << ": error message: " << e << endl;
      QDP_abort(1);
    }
    
    pop(xml_out);  // make_source

//    // Reset the seed
//    QDP::RNG::setrn(ran_seed);

    snoop.stop();
    QDPIO::cout << InlineMakeSourceFermEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineMakeSourceFermEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
