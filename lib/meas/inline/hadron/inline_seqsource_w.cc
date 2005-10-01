// $Id: inline_seqsource_w.cc,v 2.2 2005-10-01 17:01:39 edwards Exp $
/*! \file
 * \brief Inline construction of sequential sources
 *
 * Sequential source construction
 */

#include "meas/inline/hadron/inline_seqsource_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/sink_smear2.h"
#include "meas/hadron/wall_qprop_w.h"
#include "meas/hadron/seqsrc_funcmap_w.h"
#include "meas/hadron/hadseqsrc_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineSeqSourceEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSeqSource(InlineSeqSourceParams(xml_in, path));
    }

    bool registerAll()
    {
      bool foo = true;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      foo &= SeqSourceCallMapEnv::registered;
      return foo;
    }

    const std::string name = "SEQSOURCE";
    const bool registered = registerAll();
  };


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineSeqSourceParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_ids", input.prop_ids);
    read(inputtop, "seqsource_id", input.seqsource_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineSeqSourceParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "prop_ids", input.prop_ids);
    write(xml, "seqsource_id", input.seqsource_id);

    pop(xml);
  }



  // Param stuff
  InlineSeqSourceParams::InlineSeqSourceParams() { frequency = 0; }

  InlineSeqSourceParams::InlineSeqSourceParams(XMLReader& xml_in, const std::string& path) 
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
  InlineSeqSourceParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "PropSink", sink_header);
    Chroma::write(xml_out, "NamedObject", named_obj);
    
    pop(xml_out);
  }


  // Function call
  void 
  InlineSeqSource::operator()(const multi1d<LatticeColorMatrix>& u,
			      XMLBufferWriter& gauge_xml,
			      unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "seqsource");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineSeqSourceEnv::name << ": propagator sequential source constructor" << endl;
    StopWatch swatch;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

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
      QDPIO::cerr << InlineSeqSourceEnv::name << ": sanity error: " << endl;
      QDP_abort(1);
    }

    //
    // Read the quark propagator and extract headers
    //
    multi1d<LatticePropagator> forward_props(params.named_obj.prop_ids.size());
    multi1d<ChromaProp_t> prop_header(params.named_obj.prop_ids.size());
    multi1d<PropSource_t> source_header(params.named_obj.prop_ids.size());
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
	// Also pull out the id of this source
	{
	  read(prop_record_xml, "/Propagator/ForwardProp", prop_header[loop]);
	  read(prop_record_xml, "/Propagator/PropSource", source_header[loop]);
	}

	// Save prop input
	write(xml_out, "Propagator_info", prop_record_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineSeqSourceEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineSeqSourceEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      pop(xml_out);
    }
    pop(xml_out);

    QDPIO::cout << "Forward propagator successfully read and parsed" << endl;

    // Derived from input prop
    int  j_decay = source_header[0].j_decay;
    multi1d<int> t_source = source_header[0].t_source;

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
    if (params.param.t_sink < 0 || params.param.t_sink >= params.param.nrow[j_decay]) 
    {
      QDPIO::cerr << "Sink time coordinate incorrect." << endl;
      QDPIO::cerr << "t_sink = " << params.param.t_sink << endl;
      QDP_abort(1);
    }

    // Only support simple s-wave states
    if (params.sink_header.wave_state != WAVE_TYPE_S_WAVE)
    {
      QDPIO::cerr << "Only support simple s-wave states" << endl;
      QDP_abort(1);
    }



    //------------------ Start main body of calculations -----------------------------

    // Do the sink smearing BEFORE the interpolating operator
    for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
    {
      switch (params.sink_header.sink_type)
      {
      case SNK_TYPE_SHELL_SINK:
      {
	QDPIO::cout << "SeqSource: do shell sink smearing" << endl;

	sink_smear2(u, forward_props[loop], 
		    params.sink_header.sinkSmearParam.wvf_kind, 
		    params.sink_header.sinkSmearParam.wvf_param, 
		    params.sink_header.sinkSmearParam.wvfIntPar, 
		    j_decay);
      }
      break;

      case SNK_TYPE_WALL_SINK:
      {
	QDPIO::cout << "SeqSource: do wall sink smearing" << endl;

	for(int i=0; i < params.param.sink_mom.size(); ++i)
	  if (params.param.sink_mom[i] != 0)
	  {
	    QDPIO::cerr << "Do not currently support non-zero momenta wall-sink smearing" << endl;
	    QDP_abort(1);
	  }

	LatticePropagator tmp_prop = forward_props[loop];
	wall_qprop(forward_props[loop], tmp_prop, phases);
      }
      break;

      default:
	break;
      }
    }
    
  
    //
    // Construct the sequential source
    //
    QDPIO::cout << "Sequential source = " << params.param.seq_src << endl;

    swatch.reset();
    swatch.start();
    LatticePropagator quark_prop_src = 
      hadSeqSource(forward_props, 
		   params.param.t_sink, 
		   params.param.sink_mom, 
		   j_decay, 
		   params.param.seq_src);
    swatch.stop();
    
    QDPIO::cout << "Hadron sequential source computed: time= " 
		<< swatch.getTimeInSeconds() 
		<< " secs" << endl;

    // Do the sink smearing AFTER the interpolating operator
    switch (params.sink_header.sink_type)
    {
    case SNK_TYPE_SHELL_SINK:
    {
      sink_smear2(u, quark_prop_src, 
		  params.sink_header.sinkSmearParam.wvf_kind, 
		  params.sink_header.sinkSmearParam.wvf_param, 
		  params.sink_header.sinkSmearParam.wvfIntPar, 
		  j_decay);
    }
    break;

    case SNK_TYPE_WALL_SINK:
    {
//    LatticePropagator tmp_prop = quark_prop_src;
//    wall_qprop(quark_prop_src, tmp_prop, phases);

      int t_sink = params.param.t_sink;

      // Project propagator onto zero momentum: Do a slice-wise sum.
      DPropagator dprop_slice = sum(quark_prop_src, phases.getSet()[t_sink]);
      quark_prop_src[phases.getSet()[t_sink]] = dprop_slice;
    }
    break;

    default:
      break;
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
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      // Sequential source header
      // Header composed of all forward prop headers
      SequentialSource_t src;
      src.sink_header = params.sink_header;
      src.seqsource_header = params.param;
      src.forward_props.resize(params.named_obj.prop_ids.size());
      for(int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
      {
	src.forward_props[loop].sink_header = params.sink_header;
	src.forward_props[loop].prop_header = prop_header[loop];
	src.forward_props[loop].source_header = source_header[loop];
      }

      XMLBufferWriter record_xml;
      push(record_xml, "SequentialSource");
      write(record_xml, ".", src);
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);  // SequentialSource

      // Store the seqsource
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.seqsource_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.seqsource_id) = quark_prop_src;
      TheNamedObjMap::Instance().get(params.named_obj.seqsource_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.seqsource_id).setRecordXML(record_xml);

      QDPIO::cout << "Sequential source successfully stored"  << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineSeqSourceEnv::name << ": dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineSeqSourceEnv::name << ": error storing seqsource: " << e << endl;
      QDP_abort(1);
    }

    pop(xml_out);    // seqsource

    QDPIO::cout << "Seqsource ran successfully" << endl;

    END_CODE();
  } 

};
