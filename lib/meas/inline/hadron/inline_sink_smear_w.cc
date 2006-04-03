// $Id: inline_sink_smear_w.cc,v 3.0 2006-04-03 04:59:02 edwards Exp $
/*! \file
 * \brief Inline construction of sink_smear
 *
 * Sink smear propagators
 */

#include "handle.h"
#include "meas/inline/hadron/inline_sink_smear_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"

#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"

namespace Chroma 
{ 
  namespace InlineSinkSmearEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSinkSmear(InlineSinkSmearParams(xml_in, path));
    }

    bool registerAll()
    {
      bool foo = true;
      foo &= QuarkSinkSmearingEnv::registered;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      return foo;
    }

    const std::string name = "SINK_SMEAR";
    const bool registered = registerAll();
  };


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSinkSmearParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    input.gauge_id = InlineDefaultGaugeField::readGaugeId(inputtop, "gauge_id");
    read(inputtop, "prop_id", input.prop_id);
    read(inputtop, "smeared_prop_id", input.smeared_prop_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSinkSmearParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_id", input.prop_id);
    write(xml, "smeared_prop_id", input.smeared_prop_id);

    pop(xml);
  }


  // Param stuff
  InlineSinkSmearParams::InlineSinkSmearParams() { frequency = 0; }

  InlineSinkSmearParams::InlineSinkSmearParams(XMLReader& xml_in, const std::string& path) 
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
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineSinkSmearParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  void 
  InlineSinkSmear::operator()(unsigned long update_no,
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
      QDPIO::cerr << InlineSinkSmearEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSinkSmearEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "sink_smear");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineSinkSmearEnv::name << ": Sink smearing for propagators" << endl;

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //
    // Read the quark propagator and extract headers
    //
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to read forward propagator" << endl;
    try
    {
      // Grab a copy of the propagator. Will modify it later.
      quark_propagator = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	
      // Snarf the prop info. This is will throw if the prop_id is not there
      XMLReader prop_file_xml, prop_record_xml;
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_record_xml);

      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      {
	read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
	read(prop_record_xml, "/Propagator/PropSource", source_header);
      }
    }    
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineSinkSmearEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSinkSmearEnv::name << ": error extracting prop_header: " << e << endl;
      QDP_abort(1);
    }

    // Derived from input prop
    int  j_decay = source_header.j_decay;

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						   phases.getSet());

      push(xml_out, "Forward_prop_correlator");
      write(xml_out, "forward_prop_corr", forward_prop_corr);
      pop(xml_out);
    }


    //
    // Sink smear the propagator
    //
    try
    {
      QDPIO::cout << "Sink_xml = " << params.param.sink << endl;

      std::istringstream  xml_s(params.param.sink);
      XMLReader  sinktop(xml_s);
      const string sink_path = "/Sink";
      QDPIO::cout << "Sink = " << params.param.sink_type << endl;

      Handle< QuarkSourceSink<LatticePropagator> >
	sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(params.param.sink_type,
									 sinktop,
									 sink_path,
									 u));
      (*sinkSmearing)(quark_propagator);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << InlineSinkSmearEnv::name << ": Caught Exception creating sink: " << e << endl;
      QDP_abort(1);
    }


    /*
     * Save sink smeared propagator
     */
    // Save some info
    write(xml_out, "PropSource", source_header);
    write(xml_out, "ForwardProp", prop_header);
    write(xml_out, "PropSink", params.param);

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator), 
					   phases.getSet());

      push(xml_out, "SinkSmearedProp_correlator");
      write(xml_out, "sink_smeared_prop_corr", prop_corr);
      pop(xml_out);
    }


    // Save the propagator
    try
    {
      XMLBufferWriter file_xml;
      push(file_xml, "sink_smear");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      push(record_xml, "SinkSmear");
      write(record_xml, "PropSink", params.param);
      write(record_xml, "ForwardProp", prop_header);
      write(record_xml, "PropSource", source_header);
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);

      // Write the smeared prop xml info
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.smeared_prop_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.smeared_prop_id) 
	= quark_propagator;
      TheNamedObjMap::Instance().get(params.named_obj.smeared_prop_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.smeared_prop_id).setRecordXML(record_xml);

      QDPIO::cout << "Sink successfully updated" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineSinkSmearEnv::name << ": dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSinkSmearEnv::name << ": error message: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);  // sink_smear


    snoop.stop();
    QDPIO::cout << InlineSinkSmearEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineSinkSmearEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
