// $Id: seqprop_collect.cc,v 1.1 2004-04-23 15:54:05 bjoo Exp $
/*! \file
 *  \brief Main code for sequential propagator generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;


/*
 * Input 
 */
//! Parameters for running program
struct Param_t
{
  InvertParam_t    invParam;

  bool             nonRelSeqProp;
  multi1d<int>     Seq_src;    // integer array holding sequential source numbers
  multi1d<int>     sink_mom;
  int              t_sink;

  multi1d<int>     nrow;
};


//! Propagators
struct Prop_t
{
  string           prop_file;  // The files is expected to be in SciDAC format!
  multi1d<string>  seqprop_files;  // The files is expected to be in SciDAC format!
  QDP_volfmt_t     seqprop_volfmt;
};


//! Mega-structure of all input
struct Seqprop_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_file", input.prop_file);
  read(inputtop, "seqprop_files", input.seqprop_files);
  read(inputtop, "seqprop_volfmt", input.seqprop_volfmt);
}


//! Reader for input parameters
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 2:
    param.nonRelSeqProp = false;
    break;

    /**************************************************************************/
  case 3:
    read(paramtop, "nonRelSeqProp", param.nonRelSeqProp);
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "Seq_src", param.Seq_src);
  read(paramtop, "InvertParam", param.invParam);

  read(paramtop, "t_sink", param.t_sink);
  read(paramtop, "sink_mom", param.sink_mom);

  read(paramtop, "nrow", param.nrow);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Seqprop_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the forward_prop/seqprop info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}


//! Sequential propagator generation
/*
 *  \defgroup propagator Propagator generation
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

  START_CODE("seqprop");

  // Input parameter structure
  Seqprop_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/seqpropComponent", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Sanity checks
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;
  
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

  // Read in the configuration along with relevant information.
  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "seqpropCollect");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 2);
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
  XMLReader prop_file_xml, prop_record_xml;
  LatticePropagator quark_propagator;
  ChromaProp_t prop_header;
  PropSource_t source_header;
  {
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
  }

  // Derived from input prop
  int  j_decay = source_header.j_decay;
  multi1d<int> boundary = prop_header.boundary;
  multi1d<int> t_source = source_header.t_source;

  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, j_decay);
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


  // Make a sink propagator header from the source header
  // This clones params in the correct way - we use the same sink
  // smearing as we did at the source
  PropSink_t  sink_header;
  initHeader(sink_header, source_header);
  
  // Determine what kind of sink to use
  bool Sl_snk = (sink_header.sink_type == SNK_TYPE_SHELL_SINK) ? true : false;

  // Save prop input
  write(xml_out, "PropSource", source_header);
  write(xml_out, "ForwardProp", prop_header);
  write(xml_out, "PropSink", sink_header);

  // Read the quark propagator and extract headers
  //
  XMLArrayWriter  xml_seq_src(xml_out, input.param.Seq_src.size());
  push(xml_seq_src, "Sequential_source");

  for(int seq_src_ctr = 0; seq_src_ctr < input.param.Seq_src.size(); seq_src_ctr++)
  {
    QDPIO::cout << "collecting components  for seq_src number = " 
		<< seq_src_ctr << endl;

    push(xml_seq_src);
    write(xml_seq_src, "seq_src_ctr", seq_src_ctr);

    int seq_src_value = input.param.Seq_src[seq_src_ctr]; /* Assign the particular 
							     source type */

    // Allocate space for the sequential source
    LatticePropagator seq_quark_prop;
    LatticeFermion psi;
    int max_spin;
    if( input.param.nonRelSeqProp ) { 
      max_spin = Ns/2;
    }
    else {
      max_spin = Ns;
    }

    XMLReader file_xml_in, record_xml_in;
    for(int spin = 0; spin < max_spin; spin++) {
      for(int color =0; color < Nc; color++) { 
	
	ostringstream infile;
	infile <<  input.prop.seqprop_files[seq_src_ctr] << "_component_s" << spin << "_c" << color;

	readFermion(file_xml_in, record_xml_in, psi, infile.str(), 
		    QDPIO_SERIAL);

	FermToProp(psi, seq_quark_prop, color, spin);
      }
    }

    // Now trick of negatings spin components max_spin -> Ns 
    // this is the effect of a gamma_5 
    if( input.param.nonRelSeqProp ) { 
      for(int spin = max_spin; spin < Ns; spin++) {
	for(int color =0; color < Nc; color++) { 
	  
	  PropToFerm(seq_quark_prop,psi, color, spin-max_spin);
	  
	  FermToProp((-psi), seq_quark_prop, color, spin);
	}
      }
    }

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      multi1d<Double> backward_prop_corr = sumMulti(localNorm2(seq_quark_prop), 
						    phases.getSet());
	
      push(xml_seq_src, "Backward_prop_correlator");
      write(xml_seq_src, "backward_prop_corr", backward_prop_corr);
      pop(xml_seq_src);
    }


    /*
     *  Write the sequential propagator out to disk
     */
    {
      XMLBufferWriter file_xml;
      push(file_xml, "seqprop");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      // Make a seqprop structure
      ChromaSeqProp_t seqprop_header;
      seqprop_header.invParam = input.param.invParam;
      seqprop_header.nonRelSeqProp  = input.param.nonRelSeqProp;
      seqprop_header.Seq_src  = seq_src_value;
      seqprop_header.sink_mom = input.param.sink_mom;
      seqprop_header.t_sink   = input.param.t_sink;
      seqprop_header.nrow     = input.param.nrow;

      XMLBufferWriter record_xml;
      push(record_xml, "SeqProp");
      write(record_xml, "SequentialProp", seqprop_header);
      write(record_xml, "PropSink", sink_header);
      write(record_xml, "ForwardProp", prop_header);
      write(record_xml, "PropSource", source_header);
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);

      // Write the seqprop
      writeQprop(file_xml, record_xml, seq_quark_prop,
		 input.prop.seqprop_files[seq_src_ctr], 
		 input.prop.seqprop_volfmt, QDPIO_SERIAL);
    }

    /*
     *  In the case of the pion, we know that the exponentiated propagator
     *  back to the source should be the pion correlator at time-slice
     *  zero, and so will write this out
     */

    if(seq_src_value == 10)
    {
      Complex pion_src;
      seqPionTest(pion_src, seq_quark_prop, t_source);
	
      push(xml_seq_src,"Seq_propagator_test");
      write(xml_seq_src, "pion_src", pion_src);
      pop(xml_seq_src);
    }

    pop(xml_seq_src); /// Elem

  } /* end loop over sequential sources */
      
  pop(xml_seq_src);  // Sequential_source

  pop(xml_out);    // seqprop

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE("seqprop");

  exit(0);
}

