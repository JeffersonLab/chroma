// $Id: bar3ptfn.cc,v 1.31 2004-04-17 03:36:36 edwards Exp $
/*! \file
 * \brief Main program for computing 3pt functions
 *
 * Main program for computing 3pt functions
 */

#include "chroma.h"

using namespace QDP;


/*
 * Input 
 */
//! Parameters for running program
struct Param_t
{
  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  multi1d<int> nrow;
};

//! Propagators
struct Prop_t
{
  string           prop_file;  // The files is expected to be in SciDAC format!
  multi1d<string>  seqprop_files;  // The files is expected to be in SciDAC format!
};


//! Mega-structure of all input
struct Bar3ptfn_input_t
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
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 6:
    /**************************************************************************/
    break;

  default :
    /**************************************************************************/

    QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "mom2_max", param.mom2_max);
  read(paramtop, "nrow", param.nrow);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Bar3ptfn_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the propagator(s) info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading prop data: " << e << endl;
    QDP_abort(1);
  }
}


//--------------------------------------------------------------
struct Output_version_t
{
  int out_version;
};

struct FormFac_sequential_source_t
{
  int               seq_src_value;
  multi1d<int>      t_source;
  int               t_sink;
  multi1d<int>      sink_mom;
  Complex           seq_hadron_0;
  FormFac_insertions_t  formFacs;
};

struct FormFac_Wilson_3Pt_fn_measurements_t
{
  int  output_version;   // Unique id for each output version of the structures
  multi1d<FormFac_sequential_source_t> seqsrc;
};

struct Bar3ptfn_t
{
  Output_version_t  output_version;  // carry over from nml/xml like output versons
  Param_t           param;
  FormFac_Wilson_3Pt_fn_measurements_t  bar;
};


// params
void write(BinaryWriter& bin, const Output_version_t& ver)
{
  write(bin, ver.out_version);
}

// params
void write(BinaryWriter& bin, const Param_t& param)
{
  write(bin, param.mom2_max);
  write(bin, param.nrow);
}


// 
void write(BinaryWriter& bin, const FormFac_sequential_source_t& src)
{
  write(bin, src.seq_src_value);
  write(bin, src.t_source);
  write(bin, src.t_sink);
  write(bin, src.sink_mom);
  write(bin, src.seq_hadron_0);
  write(bin, src.formFacs);
}

// Write a hadron measurements
void write(BinaryWriter& bin, const FormFac_Wilson_3Pt_fn_measurements_t& had)
{
  write(bin, had.output_version);
  write(bin, had.seqsrc);
}

// Write all formfactor measurements
void write(BinaryWriter& bin, const Bar3ptfn_t& bar)
{
  write(bin, bar.output_version);
  write(bin, bar.param);
  write(bin, bar.bar);
}

//--------------------------------------------------------------





//! Main program for computing 3pt functions
/*! Main program */
int
main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  Bar3ptfn_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/bar3ptfn", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl;
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;
  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  // Read in the configuration along with relevant information.
  QDPIO::cout << "Attempt to initialize the gauge field" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;
  case CFG_TYPE_SZINQIO:
    readGauge(gauge_file_xml, gauge_xml, u, input.cfg.cfg_file, QDPIO_SERIAL);
    break;
  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDPIO::cerr << "Configuration type is unsupported." << endl;
    QDP_abort(1);
  }

  // Next check the gauge field configuration by reunitarizing.
  unitarityCheck(u);

  QDPIO::cout << "Gauge field successfully initialized" << endl;


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "bar3ptfn");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 9);
  pop(xml_out);

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
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
      QDP_abort(1);
    }
  }
  QDPIO::cout << "Forward propagator successfully read" << endl;

  // Derived from input prop
  int  j_decay = source_header.j_decay;
  multi1d<int> t_source = source_header.t_source;

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

  // Determine what kind of source to use
  bool Sl_src = (source_header.source_type == SRC_TYPE_SHELL_SOURCE) ? true : false;

  // Save prop input
  write(xml_out, "ForwardProp", prop_header);
  write(xml_out, "PropSource", source_header);

  write(xml_out,"Sl_src",Sl_src);


  // Big nested structure that is image of entire file
  Bar3ptfn_t  bar3pt;
  bar3pt.output_version.out_version = 9;  // bump this up everytime something changes
  bar3pt.param = input.param; // copy entire structure

  push(xml_out, "Wilson_3Pt_fn_measurements");

  // Big nested structure that is image of all form-factors
//    FormFac_Wilson_3Pt_fn_measurements_t  formfacs;
  bar3pt.bar.output_version = 2;  // bump this up everytime something changes
  bar3pt.bar.seqsrc.resize(input.prop.seqprop_files.size());

  XMLArrayWriter  xml_seq_src(xml_out, input.prop.seqprop_files.size());
  push(xml_seq_src, "Sequential_source");

  for (int seq_src_ctr = 0; seq_src_ctr < input.prop.seqprop_files.size(); ++seq_src_ctr) 
  {
    push(xml_seq_src);
    write(xml_seq_src, "seq_src_ctr", seq_src_ctr);

    // Read the sequential propagator
    // Read the quark propagator and extract headers
    LatticePropagator seq_quark_prop;
    ChromaSeqProp_t seqprop_header;
    PropSink_t sink_header;
    {
      XMLReader seqprop_file_xml, seqprop_record_xml;
      readQprop(seqprop_file_xml, 
		seqprop_record_xml, seq_quark_prop,
		input.prop.seqprop_files[seq_src_ctr], QDPIO_SERIAL);

      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      // NEED SECURITY HERE - need a way to cross check props. Use the ID.
      try
      {
	read(seqprop_record_xml, "/SeqProp/SequentialProp", seqprop_header);
	read(seqprop_record_xml, "/SeqProp/PropSink", sink_header);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting seqprop header: " << e << endl;
	QDP_abort(1);
      }
    }
    QDPIO::cout << "Sequential propagator successfully read" << endl;

    // Save seqprop input
    write(xml_seq_src, "SequentialProp", seqprop_header);
    write(xml_seq_src, "PropSink", sink_header);

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);
      
      multi1d<Double> backward_prop_corr = sumMulti(localNorm2(seq_quark_prop), 
						    phases.getSet());
      
      push(xml_seq_src, "Backward_prop_correlator");
      write(xml_seq_src, "backward_prop_corr", backward_prop_corr);
      pop(xml_seq_src);
    }

    // Derived from input seqprop
    int seq_src_value = seqprop_header.Seq_src;
    int           t_sink   = seqprop_header.t_sink;
    multi1d<int>  sink_mom = seqprop_header.sink_mom;

    // Output is driven by the type of sequential propagator
    if ((0 <= seq_src_value) && (seq_src_value <= 9)) {
      write(xml_seq_src, "hadron_type", "BARYON");
    } else if ((10 <= seq_src_value) && (seq_src_value <= 20)) {
      write(xml_seq_src, "hadron_type", "MESON");
    } else if ((21 <= seq_src_value) && (seq_src_value <= 30)) {
      write(xml_seq_src, "hadron_type", "BARYON");
    } else {
      QDP_error_exit("Unknown sequential source type", seq_src_value);
    }

    write(xml_seq_src, "seq_src_value", seq_src_value);
    write(xml_seq_src, "t_source", t_source);
    write(xml_seq_src, "t_sink", t_sink);
    write(xml_seq_src, "sink_mom", sink_mom);
	
    bar3pt.bar.seqsrc[seq_src_ctr].seq_src_value = seq_src_value;
    bar3pt.bar.seqsrc[seq_src_ctr].t_source      = t_source;
    bar3pt.bar.seqsrc[seq_src_ctr].t_sink        = t_sink;
    bar3pt.bar.seqsrc[seq_src_ctr].sink_mom      = sink_mom;
	
//      xml_seq_src.flush();

    // Construct the two-pt function from the source point to the sink
    // using only the seq. quark prop.
    // Take hermitian conjugate of the seq. prop, multiply on both sides
    // with gamma_5 = Gamma(G5) and take the trace
    // Use indexing to pull out precisely the source point.
    int G5 = Ns*Ns-1;

    // Contract the sequential propagator with itself
    // to form the 2-pt function at the source.
    // Do "source" smearing, if needed
    {
      LatticePropagator seq_quark_prop_tmp = seq_quark_prop;

      if (Sl_src) {
	sink_smear2(u, seq_quark_prop_tmp, 
		    source_header.sourceSmearParam.wvf_kind, 
		    source_header.sourceSmearParam.wvf_param,
		    source_header.sourceSmearParam.wvfIntPar, 
		    j_decay);
      }

      // Compute the 2pt function at the source - used in later normalizations
      Complex seq_hadron_0 =
	peekSite(trace(adj(Gamma(G5)*seq_quark_prop_tmp*Gamma(G5))), t_source);

      write(xml_seq_src, "seq_hadron_0", seq_hadron_0);
      bar3pt.bar.seqsrc[seq_src_ctr].seq_hadron_0 = seq_hadron_0;
    }

    // Now the 3pt contractions
    SftMom phases(input.param.mom2_max, sink_mom, false, j_decay);
    FormFac(bar3pt.bar.seqsrc[seq_src_ctr].formFacs, 
	    u, quark_propagator, seq_quark_prop, phases, t_source[j_decay]);

    pop(xml_seq_src);   // elem
  } // end loop over sequential sources

  pop(xml_seq_src);  // Sequential_source

//    BinaryWriter  bin_out("bar3ptfn.dat");
//    write(bin_out, bar3ptfn);
//    bin_out.close();

  pop(xml_out);  // Wilson_3Pt_fn_measurements

  // Close the namelist output file XMLDAT
  pop(xml_out);     // bar3ptfn

  BinaryWriter  bin_out("bar3ptfn.dat");
  write(bin_out, bar3pt);
  bin_out.close();

  xml_in.close();
  xml_out.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
