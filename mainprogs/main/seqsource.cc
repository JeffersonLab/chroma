// $Id: seqsource.cc,v 1.1 2004-04-24 03:26:32 edwards Exp $
/*! \file
 *  \brief Main code for sequential source construction
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;


/*
 * Input 
 */

//! Propagators
struct Prop_t
{
  string           prop_file;  // The file is expected to be in SciDAC format!
  string           seqsource_file;  // The file is expected to be in SciDAC format!
  QDP_volfmt_t     seqsource_volfmt;
};


//! Mega-structure of all input
struct SeqSource_input_t
{
  SeqSource_t       param;
  PropSink_t        sink_header;
  Cfg_t             cfg;
  Prop_t            prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "seqsource_file", input.seqsource_file);
  read(inputtop, "seqsource_volfmt", input.seqsource_volfmt);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, SeqSource_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // The parameters for smearing the sink
    read(inputtop, "PropSink", input.sink_header);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the forward_prop/seqsource info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}


//! Sequential source generation
/*
 *  \defgroup seqsource Sequential source generation
 *  \ingroup main
 *
 *  Read quark propagators, compute the sequential sources needed
 *  for baryon and meson form factors and/or structure function moments.
 *
 *  This routine does not compute the form factors or moments --
 *  that is done in separate routines....
 *
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE("seqsource");

  // Input parameter structure
  SeqSource_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/seqsource", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Sanity checks
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  QDPIO::cout << "     Computing sequential source of type "
	      << input.param.seq_src << endl;
  
  QDPIO::cout << "     Volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;


  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "seqsource");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();


  //
  // Read the quark propagator and extract headers
  //
  LatticePropagator quark_propagator;
  ChromaProp_t prop_header;
  PropSource_t source_header;
  {
    XMLReader prop_file_xml, prop_record_xml;

    QDPIO::cout << "Attempt to read forward propagator" << endl;
    readQprop(prop_file_xml, 
	      prop_record_xml, quark_propagator,
	      input.prop.prop_file, QDPIO_SERIAL);
    QDPIO::cout << "Forward propagator successfully read" << endl;
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
      read(prop_record_xml, "/Propagator/PropSource", source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
      throw;
    }

    // Save prop input
    write(xml_out, "Forward_propagator", prop_record_xml);
  }

  // Derived from input prop
  int  j_decay = source_header.j_decay;
  multi1d<int> boundary = prop_header.boundary;
  multi1d<int> t_source = source_header.t_source;

  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, j_decay);

  // Sanity check - write out the norm2 of the forward prop in the j_decay direction
  // Use this for any possible verification
  {
    multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						 phases.getSet());

    push(xml_out, "Forward_prop_correlator");
    write(xml_out, "forward_prop_corr", forward_prop_corr);
    pop(xml_out);
  }

  // A sanity check
  if (input.param.t_sink < 0 || input.param.t_sink >= input.param.nrow[j_decay]) 
  {
    QDPIO::cerr << "Sink time coordinate incorrect." << endl;
    QDPIO::cerr << "t_sink = " << input.param.t_sink << endl;
    QDP_abort(1);
  }

  // Only support simple s-wave states
  if (input.sink_header.wave_state != WAVE_TYPE_S_WAVE)
  {
    QDPIO::cerr << "Only support simple s-wave states" << endl;
    QDP_abort(1);
  }



  //------------------ Start main body of calculations -----------------------------

  // Do the sink smearing BEFORE the interpolating operator
  switch (input.sink_header.sink_type)
  {
  case SNK_TYPE_SHELL_SINK:
  {
    QDPIO::cout << "SeqSource: do shell sink smearing" << endl;

    sink_smear2(u, quark_propagator, 
		source_header.sourceSmearParam.wvf_kind, 
		source_header.sourceSmearParam.wvf_param, 
		source_header.sourceSmearParam.wvfIntPar, 
		j_decay);
  }
  break;

  case SNK_TYPE_WALL_SINK:
  {
    QDPIO::cout << "SeqSource: do wall sink smearing" << endl;

    for(int i=0; i < input.param.sink_mom.size(); ++i)
      if (input.param.sink_mom[i] != 0)
      {
	QDPIO::cerr << "Do not currently support non-zero momenta wall-sink smearing" << endl;
	QDP_abort(1);
      }

    LatticePropagator tmp_prop = quark_propagator;
    wall_qprop(quark_propagator, tmp_prop, phases);
  }
  break;

  default:
    break;
  }
    
  
  //
  // Construct the sequential source
  //
  // Allocate space for the sequential source
  LatticePropagator quark_prop_src;

  /*
   *  Sources 0 -> 9 corresponding to Baryon sequential sources
   *  Sources 10 -> 19 corresponds to a Meson sequential source
   *  Souces  21 -> 29 are additional Baryon ones we thought of
   *
   *  Note that not all the source values are necessarily implemented
   *
   */

  SeqSourceType seq_src = input.param.seq_src;


  if(((0 <= seq_src) && (seq_src <= 9)) ||
     ((21 <= seq_src) && (seq_src <= 29))) 
  {
    // Computation of the Baryon sequential source
    barSeqSource(quark_propagator, quark_propagator, quark_prop_src, 
		 input.param.t_sink, 
		 input.param.sink_mom, 
		 j_decay, 
		 seq_src);
  }
  else if ((10 <= seq_src) && (seq_src <= 20))
  {
    // Computation of the Meson sequential source
    mesonSeqSource(quark_propagator, quark_prop_src, 
		   input.param.t_sink, 
		   input.param.sink_mom, 
		   j_decay, 
		   seq_src);
  }
  else{
    QDP_error_exit("Unknown sequential source type", seq_src);
  }


  // Do the sink smearing AFTER the interpolating operator
  switch (input.sink_header.sink_type)
  {
  case SNK_TYPE_SHELL_SINK:
  {
    sink_smear2(u, quark_prop_src, 
		source_header.sourceSmearParam.wvf_kind, 
		source_header.sourceSmearParam.wvf_param, 
		source_header.sourceSmearParam.wvfIntPar, 
		j_decay);
  }
  break;

  case SNK_TYPE_WALL_SINK:
  {
    LatticePropagator tmp_prop = quark_prop_src;
    wall_qprop(quark_prop_src, tmp_prop, phases);
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
   *  Write the sequential source out to disk
   */
  {
    XMLBufferWriter file_xml;
    push(file_xml, "seqsource");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(file_xml, "id", id);
    pop(file_xml);

    XMLBufferWriter record_xml;
    push(record_xml, "SequentialSource");
    write(record_xml, "SeqSourceSinkSmear", input.sink_header);
    write(record_xml, "SeqSource", input.param);

    // For now, only have 1 set of forward props to tie
    XMLArrayWriter  xml_props(record_xml, 1);
    push(xml_props, "ForwardProps");

    for (int i=0; i < 1; ++i)
    {
      push(xml_props);
      write(xml_props, "PropSink", input.sink_header);
      write(xml_props, "ForwardProp", prop_header);
      write(xml_props, "PropSource", source_header);
      pop(xml_props);
    }
    pop(xml_props);  // ForwardProps

    write(record_xml, "Config_info", gauge_xml);

    // Write the seqsource
    writeQprop(file_xml, record_xml, quark_prop_src,
	       input.prop.seqsource_file, 
	       input.prop.seqsource_volfmt, QDPIO_SERIAL);
  }

  pop(xml_out);    // seqsource

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE("seqsource");

  exit(0);
}

