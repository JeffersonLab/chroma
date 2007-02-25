// $Id: inline_seqsource_w.cc,v 3.6 2007-02-25 22:39:28 edwards Exp $
/*! \file
 * \brief Inline construction of sequential sources
 *
 * Sequential source construction
 */

#include "handle.h"
#include "meas/inline/hadron/inline_seqsource_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/hadron/seqsource_aggregate_w.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSeqSourceEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_ids", input.prop_ids);
    read(inputtop, "seqsource_id", input.seqsource_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSeqSourceEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_ids", input.prop_ids);
    write(xml, "seqsource_id", input.seqsource_id);

    pop(xml);
  }


  namespace InlineSeqSourceEnv 
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

    const std::string name = "SEQSOURCE";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {	
	success &= QuarkSinkSmearingEnv::registerAll();
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

	// The parameters holds the version number
	read(paramtop, "Param", param);

	// The parameters for smearing the sink
	read(paramtop, "PropSink", sink_header);

	// Read in the forward_prop/seqsource info
	read(paramtop, "NamedObject", named_obj);
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
      write(xml_out, "PropSink", sink_header);
      write(xml_out, "NamedObject", named_obj);
    
      pop(xml_out);
    }


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
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
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "seqsource");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator sequential source constructor" << endl;
      StopWatch swatch;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      // Sanity check
      if (params.named_obj.prop_ids.size() == 0)
      {
	QDPIO::cerr << name << ": sanity error: " << endl;
	QDP_abort(1);
      }

      //
      // Read the quark propagator and extract headers
      //
      multi1d<LatticePropagator> forward_props(params.named_obj.prop_ids.size());
      multi1d<ForwardProp_t> forward_headers(params.named_obj.prop_ids.size());
      push(xml_out, "Forward_prop_infos");
      for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
      {
	push(xml_out, "elem");
	try
	{
	  // Snarf the data into a copy
	  forward_props[loop] =
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[loop]);
	
	  // Snarf the source info. This is will throw if the source_id is not there
	  XMLReader prop_file_xml, prop_record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getFileXML(prop_file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getRecordXML(prop_record_xml);
   
	  // Try to invert this record XML into a ChromaProp struct
	  {
	    Propagator_t  header;
	    read(prop_record_xml, "/Propagator", header);

	    forward_headers[loop].prop_header   = header.prop_header;
	    forward_headers[loop].source_header = header.source_header;
	    forward_headers[loop].gauge_header  = header.gauge_header;
	  }

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
	  QDPIO::cerr << name << ": map call failed: " << e 
		      << endl;
	  QDP_abort(1);
	}
	pop(xml_out);
      }
      pop(xml_out);

      QDPIO::cout << "Forward propagator successfully read and parsed" << endl;

      // Derived from input prop
      int j_decay  = forward_headers[0].source_header.j_decay;

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      push(xml_out, "Forward_prop_correlators");
      for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
      {
	multi1d<Double> forward_prop_corr = sumMulti(localNorm2(forward_props[loop]),
						     phases.getSet());

	push(xml_out, "elem");
	write(xml_out, "forward_prop_corr", forward_prop_corr);
	pop(xml_out);
      }
      pop(xml_out);

      // A sanity check
      if (params.param.t_sink < 0 || params.param.t_sink >= QDP::Layout::lattSize()[j_decay]) 
      {
	QDPIO::cerr << "Sink time coordinate incorrect." << endl;
	QDPIO::cerr << "t_sink = " << params.param.t_sink << endl;
	QDP_abort(1);
      }


      //------------------ Start main body of calculations -----------------------------

      LatticePropagator quark_prop_src;

      try
      {
	// Sink smear the forward propagators
	// NOTE: The smearing construction is pulled outside the loop
	// for efficiency. However, I'm anticipating that we will 
	// have different smearings at the sink of the forward props.
	// In that case, the loop needs to be in inverted.
	std::istringstream  xml_s(params.sink_header.sink.xml);
	XMLReader  sinktop(xml_s);
	QDPIO::cout << "Sink = " << params.sink_header.sink.id << endl;
	
	Handle< QuarkSourceSink<LatticePropagator> >
	  sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(params.sink_header.sink.id,
									   sinktop,
									   params.sink_header.sink.path,
									   u));

	// Do the sink smearing BEFORE the interpolating operator
	for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
	{
	  forward_headers[loop].sink_header = params.sink_header;
	  (*sinkSmearing)(forward_props[loop]);
	}
    
  
	//
	// Construct the sequential source
	//
	QDPIO::cout << "Sequential source = " << params.param.seqsrc.xml << endl;

	std::istringstream  xml_seq(params.param.seqsrc.xml);
	XMLReader  seqsrctop(xml_seq);
	QDPIO::cout << "SeqSource = " << params.param.seqsrc.id << endl;
	
	Handle< HadronSeqSource<LatticePropagator> >
	  hadSeqSource(TheWilsonHadronSeqSourceFactory::Instance().createObject(params.param.seqsrc.id,
										seqsrctop,
										params.param.seqsrc.path));

	swatch.reset();
	swatch.start();
	quark_prop_src = (*hadSeqSource)(u, forward_headers, forward_props);

	swatch.stop();
    
	QDPIO::cout << "Hadron sequential source computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;


	// Do the sink smearing AFTER the interpolating operator
	(*sinkSmearing)(quark_prop_src);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in sink: " << e << endl;
	QDP_abort(1);
      }
    
    
      // Sanity check - write out the norm2 of the propagator source in the j_decay direction
      // Use this for any possible verification
      {
	multi1d<Double> seqsource_corr = sumMulti(localNorm2(quark_prop_src), 
						  phases.getSet());
	
	push(xml_out, "SeqSource_correlator");
	write(xml_out, "seqsource_corr", seqsource_corr);
	pop(xml_out);
      }


      /*
       *  Write the sequential source out to a named buffer
       */
      try
      {
	QDPIO::cout << "Attempt to store sequential source" << endl;

	XMLBufferWriter file_xml;
	push(file_xml, "seqsource");
	write(file_xml, "id", uniqueId());  // NOTE: new ID form
	pop(file_xml);

	// Sequential source header
	// Header composed of all forward prop headers
	SequentialSource_t new_header;
	new_header.sink_header      = params.sink_header;
	new_header.seqsource_header = params.param;
	new_header.forward_props    = forward_headers;
	new_header.gauge_header     = gauge_xml.printCurrentContext();

	XMLBufferWriter record_xml;
	write(record_xml, "SequentialSource", new_header);

	// Store the seqsource
	TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.seqsource_id);
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.seqsource_id) = quark_prop_src;
	TheNamedObjMap::Instance().get(params.named_obj.seqsource_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.seqsource_id).setRecordXML(record_xml);

	QDPIO::cout << "Sequential source successfully stored"  << endl;
      }
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error storing seqsource: " << e << endl;
	QDP_abort(1);
      }

      pop(xml_out);    // seqsource

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
