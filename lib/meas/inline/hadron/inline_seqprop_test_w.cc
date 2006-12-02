// $Id: inline_seqprop_test_w.cc,v 3.3 2006-12-02 18:18:07 edwards Exp $
/*! \file
 * \brief Test sequential propagator
 *
 * Sequential source test
 */

#include "handle.h"
#include "meas/inline/hadron/inline_seqprop_test_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/hadron/seqsource_aggregate_w.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sources/source_smearing_aggregate.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/make_xml_file.h"

namespace Chroma 
{ 
  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSeqPropTestEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "sink_ids", input.sink_ids);
    read(inputtop, "seqprop_id", input.seqprop_id);
    read(inputtop, "gamma_insertion", input.gamma_insertion);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSeqPropTestEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "sink_ids", input.sink_ids);
    write(xml, "seqprop_id", input.seqprop_id);
    write(xml, "gamma_insertion", input.gamma_insertion);

    pop(xml);
  }


  namespace InlineSeqPropTestEnv 
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

    const std::string name = "SEQPROP_TEST";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {	
	success &= QuarkSourceSmearingEnv::registerAll();
	success &= HadronSeqSourceEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }



    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// The parameters for smearing the sink
	read(paramtop, "PropSourceSmear", smear_header);

	// Read in the forward_prop/seqsource info
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
    
      write(xml_out, "PropSourceSmear", smear_header);
      write(xml_out, "NamedObject", named_obj);
    
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

	push(xml_out, "SeqPropTest");
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


    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "SeqPropTest");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": sequential propagator sequential test" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

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
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting gauge field: " << e 
		    << endl;
	QDP_abort(1);
      }

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      // Sanity check
      if (params.named_obj.sink_ids.size() == 0)
      {
	QDPIO::cerr << name << ": sanity error: " << endl;
	QDP_abort(1);
      }


      //
      // Read the forward propagator and extract headers
      //
      multi1d<LatticePropagator> forward_props(params.named_obj.sink_ids.size());
      multi1d<ForwardProp_t> forward_headers(params.named_obj.sink_ids.size());

      push(xml_out, "Forward_prop_infos");
      for(int loop=0; loop < params.named_obj.sink_ids.size(); ++loop)
      {
	push(xml_out, "elem");
	try
	{
	  // Snarf the data into a copy
	  forward_props[loop] =
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.sink_ids[loop]);
	
	  // Snarf the source info. This is will throw if the source_id is not there
	  XMLReader prop_file_xml, prop_record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.sink_ids[loop]).getFileXML(prop_file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.sink_ids[loop]).getRecordXML(prop_record_xml);
   
	  // Try to invert this record XML into a ChromaProp struct
	  // Also pull out the id of this source
	  read(prop_record_xml, "/SinkSmear", forward_headers[loop]);

	  // Save prop input
	  write(xml_out, "Propagator_info", prop_record_xml);
	}
	catch( std::bad_cast ) 
	{
	  QDPIO::cerr << name << ": caught dynamic cast error" 
		      << endl;
	  QDP_abort(1);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << name << ": error extracting forward props: " << e 
		      << endl;
	  QDP_abort(1);
	}
	pop(xml_out);
      }
      pop(xml_out);

      QDPIO::cout << "Forward propagators successfully read and parsed" << endl;


      //
      // Read the quark sequential propagator and extract headers
      //
      LatticePropagator seqprop;
      SequentialProp_t seqprop_header;

      try
      {
	// Snarf the data into a copy
	seqprop = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.seqprop_id);
	
	// Snarf the source info. This is will throw if the source_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(params.named_obj.seqprop_id).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.seqprop_id).getRecordXML(prop_record_xml);
   
	// Save prop input
	write(xml_out, "Seqprop_info", prop_record_xml);

	// Try to invert this record XML into a SequentialSource_t struct
	// Also pull out the id of this source
	read(prop_record_xml, "/SequentialProp", seqprop_header);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting seqprop: " << e 
		    << endl;
	QDP_abort(1);
      }

      QDPIO::cout << "Sequential propagator successfully read and parsed" << endl;

      // Derived from input prop
      int j_decay  = seqprop_header.seqsource_header.j_decay;

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      push(xml_out, "Forward_prop_correlators");
      for(int loop=0; loop < params.named_obj.sink_ids.size(); ++loop)
      {
	multi1d<Double> forward_prop_corr = sumMulti(localNorm2(forward_props[loop]),
						     phases.getSet());

	push(xml_out, "elem");
	write(xml_out, "forward_prop_corr", forward_prop_corr);
	pop(xml_out);
      }
      pop(xml_out);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	multi1d<Double> forward_prop_corr = sumMulti(localNorm2(seqprop),
						     phases.getSet());
      
	push(xml_out, "Seqprop_correlator");
	write(xml_out, "seqprop_corr", forward_prop_corr);
	pop(xml_out);
      }

      //------------------ Start main body of calculations -----------------------------

      try
      {
	// Source smear the forward propagators
	std::istringstream  xml_s(params.smear_header.source.xml);
	XMLReader  sourcetop(xml_s);
	QDPIO::cout << "Source smearing = " << params.smear_header.source.id << endl;
	
	Handle< QuarkSourceSink<LatticePropagator> >
	  sourceSmearing(ThePropSourceSmearingFactory::Instance().createObject(params.smear_header.source.id,
									       sourcetop,
									       params.smear_header.source.path,
									       u));

	// Source smear the sequential propagator
	(*sourceSmearing)(seqprop);


	//
	// Construct the sequential source constructor
	//
	QDPIO::cout << "Sequential source = " << seqprop_header.seqsource_header.seqsrc.xml << endl;

	std::istringstream  xml_seq(seqprop_header.seqsource_header.seqsrc.xml);
	XMLReader  seqsrctop(xml_seq);
	QDPIO::cout << "SeqSource = " << seqprop_header.seqsource_header.seqsrc.id << endl;
	
	Handle< HadronSeqSource<LatticePropagator> >
	  hadSeqSource(TheWilsonHadronSeqSourceFactory::Instance().createObject(seqprop_header.seqsource_header.seqsrc.id,
										seqsrctop,
										seqprop_header.seqsource_header.seqsrc.path));

	// Evaluate the seqprop back at the source - this is the 2-pt at the sink
	int gamma_insertion = params.named_obj.gamma_insertion;

	Complex source_value = hadSeqSource->tieBack(u, seqprop_header, seqprop, gamma_insertion);

	// Compute the 2-pt 
	Complex twopt_sink   = hadSeqSource->twoPtSink(u, forward_headers, forward_props, gamma_insertion);

	// Dump the results
	push(xml_out, "LoopBackTest");
	write(xml_out, "seq_src", seqprop_header.seqsource_header.seqsrc.id);
	write(xml_out, "gamma_insertion", gamma_insertion);
	write(xml_out, "j_decay", seqprop_header.forward_props[0].source_header.j_decay);
	write(xml_out, "t_srce", seqprop_header.forward_props[0].source_header.getTSrce());
	write(xml_out, "t_source", seqprop_header.forward_props[0].source_header.t_source);
	write(xml_out, "t_sink", seqprop_header.seqsource_header.t_sink);
	write(xml_out, "sink_mom", seqprop_header.seqsource_header.sink_mom);
	write(xml_out, "source_value", source_value);
	write(xml_out, "twopt_sink", twopt_sink);
	write(xml_out, "diff", Complex(source_value - twopt_sink));
	write(xml_out, "norm2_source_value", norm2(source_value));
	write(xml_out, "norm2_twopt_sink", norm2(twopt_sink));
	write(xml_out, "norm2_diff", norm2(source_value - twopt_sink));
	pop(xml_out);
    
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in sink: " << e << endl;
	QDP_abort(1);
      }
    

      pop(xml_out);    // seqprop_test

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
