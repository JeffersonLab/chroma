// $Id: dwf_seqprop.cc,v 1.2 2004-04-27 21:29:32 edwards Exp $
/*! \file
 *  \brief Main code for sequential propagator generation for domain wall
 * fermions (array variant). Should be eliminated when FermAct becomes
 * better.
 * 
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
  multi1d<SeqSourceType>  seq_src;    // integer array holding sequential source numbers
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

  read(paramtop, "seq_src", param.seq_src);
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
  read(xml_in, "/seqprop", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Sanity checks
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  for(int seq_src_ctr = 0; seq_src_ctr < input.param.seq_src.size(); seq_src_ctr++)
    QDPIO::cout << "     Computing sequential source of type "
		<< input.param.seq_src[seq_src_ctr] << endl;
  
  QDPIO::cout << "     Volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;


  // Read in the configuration along with relevant information.
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
    QDP_error_exit("Configuration type is unsupported.");
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "seqprop");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

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


  //------------------ Start main body of calculations -----------------------------
  /*
   * Construct fermionic BC. Need one for LatticeFermion and multi1d<LatticeFermion>
   * Note, the handle is on an ABSTRACT type
   */

  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));

  //
  // Initialize fermion action
  //
  WilsonTypeFermAct< multi1d<LatticeFermion> >* S_f_a_ptr = 0;

  switch (prop_header.FermActHandle->getFermActType())
  {
  case FERM_ACT_DWF:
    {
      QDPIO::cout << "FERM_ACT_DWF" << endl;
      const DWFFermActParams& dwf = 
	dynamic_cast<const DWFFermActParams&>(*(prop_header.FermActHandle)) ;
      
      S_f_a_ptr= new EvenOddPrecDWFermActArray(fbc_a,
					       dwf.chiralParam.OverMass,
					       dwf.Mass, 
					       dwf.chiralParam.N5);
    }
    break;
    
  case FERM_ACT_UNPRECONDITIONED_DWF:
    {
      QDPIO::cout << "FERM_ACT_UNPRECONDITONED_DWF" << endl;
      const DWFFermActParams& dwf = 
	dynamic_cast<const DWFFermActParams&>(*(prop_header.FermActHandle)) ;
      S_f_a_ptr = new UnprecDWFermActArray(fbc_a,
					   dwf.chiralParam.OverMass, 
					   dwf.Mass, 
					   dwf.chiralParam.N5);
    }
    break;
    
  default:
    QDPIO::cerr << "Unsupported fermion action" << endl;
    QDP_abort(1);
  }
  
  // Create a useable handle on the action
  // The handle now owns the pointer
  Handle< WilsonTypeFermAct< multi1d<LatticeFermion> > > S_f_a(S_f_a_ptr);

  QDPIO::cout << "Seqprop: fermion action initialized" << endl;


  if (Sl_snk);
  {
    QDPIO::cout << "Seqprop: do sink smearing" << endl;
    
    // Do the sink smearing BEFORE the interpolating operator
    sink_smear2(u, quark_propagator, 
		source_header.sourceSmearParam.wvf_kind, 
		source_header.sourceSmearParam.wvf_param, 
		source_header.sourceSmearParam.wvfIntPar, 
		j_decay);
  }
  
  //
  // Loop over the sequential propagators
  //
  XMLArrayWriter  xml_seq_src(xml_out, input.param.seq_src.size());
  push(xml_seq_src, "Sequential_source");
  
  int ncg_had = 0;			// Initialise iteration counter
  for(int seq_src_ctr = 0; seq_src_ctr < input.param.seq_src.size(); seq_src_ctr++)
  {
    push(xml_seq_src);
    write(xml_seq_src, "seq_src_ctr", seq_src_ctr);

    QDPIO::cout << "Start seqprop calculation for seq_src number = " 
		<< seq_src_ctr << endl;

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

    int seq_src_value = input.param.seq_src[seq_src_ctr]; /* Assign the particular 
							     source type */


    if(((0 <= seq_src_value) && (seq_src_value <= 9)) ||
       ((21 <= seq_src_value) && (seq_src_value <= 29))) 
    {
      // Computation of the Baryon sequential source
      barSeqSource(quark_propagator, quark_propagator, quark_prop_src, 
		   input.param.t_sink, 
		   input.param.sink_mom, 
		   j_decay, 
		   seq_src_value);
    }
    else if ((10 <= seq_src_value) && (seq_src_value <= 20))
    {
      // Computation of the Meson sequential source
      mesonSeqSource(quark_propagator, quark_prop_src, 
		     input.param.t_sink, 
		     input.param.sink_mom, 
		     j_decay, 
		     seq_src_value);
    }
    else{
      QDP_error_exit("Unknown sequential source type", seq_src_value);
    }

    if (Sl_snk)
    {
      // Do the sink smearing AFTER the interpolating operator
      sink_smear2(u, quark_prop_src, 
		  source_header.sourceSmearParam.wvf_kind, 
		  source_header.sourceSmearParam.wvf_param, 
		  source_header.sourceSmearParam.wvfIntPar, 
		  j_decay);
    }

    /*
     *  Compute the full propagator.
     */
    LatticePropagator seq_quark_prop = zero;
    {
      Handle<const ConnectState> state(S_f_a->createState(u));// inserts any BC
      int n_count;

      quarkProp4(seq_quark_prop, xml_seq_src, quark_prop_src,
		 *S_f_a, state, 
		 input.param.invParam.invType, 
		 input.param.invParam.RsdCG, 
		 input.param.invParam.MaxCG,
		 input.param.nonRelSeqProp,
		 n_count);

      ncg_had += n_count;
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
      seqprop_header.seq_src  = input.param.seq_src[seq_src_ctr];
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

    pop(xml_seq_src);   // elem
  } /* end loop over sequential sources */
      
  pop(xml_seq_src);  // Sequential_source

  push(xml_out,"Relaxation_Iterations");
  write(xml_out, "ncg_had", ncg_had);
  pop(xml_out);

  pop(xml_out);    // seqprop

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE("seqprop");

  exit(0);
}

