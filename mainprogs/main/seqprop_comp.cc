// $Id: seqprop_comp.cc,v 1.4 2004-04-28 14:56:11 bjoo Exp $
/*! \file
 *  \brief Main code for sequential propagator generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;

struct Component_t { 
  int color;
  int spin;
};

/*
 * Input 
 */
//! Parameters for running program
struct Param_t
{
  InvertParam_t    invParam;

  bool             nonRelSeqProp;
  multi1d<SeqSourceType>     seq_src;    // integer array holding sequential source numbers
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
struct SeqpropComponent_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
  multi1d<Component_t> components;
};


void read(XMLReader& xml, const string& path, Component_t &comp)
{
  XMLReader top(xml,path);
  try {
    read(top, "color", comp.color);
    read(top, "spin",  comp.spin);
  }
  catch( const string& e ) {
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }  
  if( comp.color < 0 || comp.color >= Nc ) { 
    QDPIO::cerr << "Component color >= Nc. color = " << comp.color << endl;
    QDP_abort(1);
  }

  if( comp.spin < 0 || comp.spin >= Ns ) { 
    QDPIO::cerr << "Component spin >= Ns.  spin = " << comp.spin << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const Component_t &comp)
{
  
  push( xml, path );

  write(xml, "color", comp.color);
  write(xml, "spin",  comp.spin);
  
  pop( xml );
}


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

  int Seq_src;
  switch (version) 
  {
    /**************************************************************************/
  case 2:
    param.nonRelSeqProp = false;
    read(paramtop, "Seq_src", Seq_src);
    param.seq_src = SeqSourceType(Seq_src);
    break;

    /**************************************************************************/
  case 3:
    read(paramtop, "nonRelSeqProp", param.nonRelSeqProp);
    read(paramtop, "Seq_src", Seq_src);
    param.seq_src = SeqSourceType(Seq_src);
    break;

    /**************************************************************************/
  case 4:
    read(paramtop, "nonRelSeqProp", param.nonRelSeqProp);
    read(paramtop, "seq_src", param.seq_src);
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "InvertParam", param.invParam);

  read(paramtop, "t_sink", param.t_sink);
  read(paramtop, "sink_mom", param.sink_mom);

  read(paramtop, "nrow", param.nrow);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, SeqpropComponent_input_t& input)
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

     read(inputtop, "Components", input.components);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }

  // Sanity check. If non rel flags are given, spin components must
  // 
  for(int comp=0; comp < input.components.size(); comp++) { 
    if( input.param.nonRelSeqProp && (input.components[comp].spin >= (Ns/2))) {
      QDPIO::cerr << "nonRelSeqProp specified, but also spin component " << input.components[comp].spin << endl;
      QDPIO::cerr << "if nonRelSeqProp is specified spins must be < " << Ns/2<< endl;
      QDP_abort(1);
    }
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
  SeqpropComponent_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/seqpropComponent", input);

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

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);



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
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(boundary));
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));

  //
  // Initialize fermion action
  //
  FermionAction<LatticeFermion>* S_f_ptr = 0;
  FermionAction< multi1d<LatticeFermion> >* S_f_a_ptr = 0;

  switch (prop_header.FermActHandle->getFermActType() )
  {
  case FERM_ACT_WILSON:
  {
    const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams&>(*(prop_header.FermActHandle));

    QDPIO::cout << "FERM_ACT_WILSON" << endl;
    S_f_ptr = new EvenOddPrecWilsonFermAct(fbc, wils.Mass,
					   wils.anisoParam);
  }
  break;

  case FERM_ACT_UNPRECONDITIONED_WILSON:
  {
    const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams&>(*(prop_header.FermActHandle));

    QDPIO::cout << "FERM_ACT_UNPRECONDITIONED_WILSON" << endl;
    S_f_ptr = new UnprecWilsonFermAct(fbc, wils.Mass);
  }
  break;

  case FERM_ACT_ZOLOTAREV_4D:
  {
    QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
    const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(prop_header.FermActHandle));
    
    // Construct Fermact -- now uses constructor from the zolo4d params
    // struct
    S_f_ptr = new Zolotarev4DFermAct(fbc, zolo4d, xml_out);
  }
  break;

  case FERM_ACT_ZOLOTAREV_5D:
  {
    QDPIO::cout << "FERM_ACT_ZOLOTAREV_5D" << endl;
    const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(prop_header.FermActHandle));
    
    // Construct Fermact -- now uses constructor from the zolo4d params
    // struct
    S_f_a_ptr = new Zolotarev5DFermActArray(fbc_a, fbc, zolo5d, xml_out);
  }
  break;

  case FERM_ACT_DWF:
  {
    const DWFFermActParams& dwf = dynamic_cast<const DWFFermActParams&>(*(prop_header.FermActHandle));

    QDPIO::cout << "FERM_ACT_DWF" << endl;
    S_f_a_ptr = new EvenOddPrecDWFermActArray(fbc_a,
					      dwf.chiralParam.OverMass, 
					      dwf.Mass, 
					      dwf.chiralParam.N5);
  }
  break;

  case FERM_ACT_UNPRECONDITIONED_DWF:
  {
    const DWFFermActParams& dwf = dynamic_cast<const DWFFermActParams&>(*(prop_header.FermActHandle));

    QDPIO::cout << "FERM_ACT_UNPRECONDITONED_DWF" << endl;
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
  Handle< FermionAction<LatticeFermion> > S_f(S_f_ptr);
  Handle< FermionAction< multi1d<LatticeFermion> > > S_f_a(S_f_a_ptr);

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

    // This is a little non trivial

    // FIrst we have to set up the state -- this is fermact dependent
    const ConnectState *state_ptr;


    switch(prop_header.FermActHandle->getFermActType()) {
    case FERM_ACT_WILSON:
      state_ptr = S_f->createState(u);
      break;
    case FERM_ACT_UNPRECONDITIONED_WILSON:
      state_ptr = S_f->createState(u);
      break;
    case FERM_ACT_DWF:
      state_ptr = S_f->createState(u);
      break;
    case FERM_ACT_UNPRECONDITIONED_DWF:
      state_ptr = S_f->createState(u);
      break;

    case FERM_ACT_ZOLOTAREV_4D:
      {    
	const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(prop_header.FermActHandle));
	const Zolotarev4DFermAct& S_zolo4 = dynamic_cast<const Zolotarev4DFermAct&>(*S_f);

	state_ptr = S_zolo4.createState(u, zolo4d.StateInfo, xml_out,zolo4d.AuxFermActHandle->getMass());
       
      }
      break;
    case FERM_ACT_ZOLOTAREV_5D:
      {
	const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(prop_header.FermActHandle));
	const Zolotarev5DFermActArray& S_zolo5 = dynamic_cast<const Zolotarev5DFermActArray&>(*S_f_a);


	state_ptr = S_zolo5.createState(u, zolo5d.StateInfo, xml_out, zolo5d.AuxFermActHandle->getMass());
      }
      break;

    default:
      QDPIO::cerr << "Unsupported fermion action (state creation)" << endl;
      QDP_abort(1);
    }
      
    // Now do the quarkprop here again depending on which of the
    // action pointers is not null
    
    Handle<const ConnectState> state(state_ptr);  // inserts any BC
    
    LatticeFermion seq_prop_component = zero;
    LatticeFermion seq_src_component = zero;

    for(int comp = 0; comp < input.components.size(); comp++) { 

      PropToFerm(quark_prop_src, seq_src_component, input.components[comp].color,
		 input.components[comp].spin);

      seq_prop_component = zero;

      Real fact = Real(1) / sqrt(norm2(seq_src_component));
      seq_src_component *= fact;
      
      int n_count = 0;
      
      if( S_f_ptr != 0x0 ) { 
	
	S_f->qprop(seq_prop_component,
		   state,
		   seq_src_component,
		   input.param.invParam.invType, 
		   input.param.invParam.RsdCG, 
		   input.param.invParam.MaxCG,
		   n_count);
	
	/*	  quarkProp4(seq_quark_prop, xml_seq_src, quark_prop_src,
	 *S_f, state, 
	 input.param.invParam.invType, 
	 input.param.invParam.RsdCG, 
	 input.param.invParam.MaxCG,
	 input.param.nonRelSeqProp,
	 n_count);
	*/
      }
      else if ( S_f_a_ptr != 0x0 ) { 
	
	S_f_a->qprop(seq_prop_component,
		     state,
		     seq_src_component,
		     input.param.invParam.invType, 
		     input.param.invParam.RsdCG, 
		     input.param.invParam.MaxCG,
		     n_count);
	
	/* quarkProp4(seq_quark_prop, xml_seq_src, quark_prop_src,
	 *S_f_a, state, 
	 input.param.invParam.invType, 
	 input.param.invParam.RsdCG, 
	 input.param.invParam.MaxCG,
	 input.param.nonRelSeqProp,
	 n_count);
	*/
      }
      else {
	QDPIO::cerr << "Both S_f_ptr and S_f_a_ptr == 0 " << endl;
	QDP_abort(1);
      }
      
      ncg_had += n_count;
      
      fact = Real(1) / fact;
      seq_prop_component *= fact;
      
      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", n_count);
      pop(xml_out);
      
      

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      
      multi1d<Double> backward_prop_corr = sumMulti(localNorm2(seq_prop_component), 
						      phases.getSet());
	
      push(xml_seq_src, "Backward_prop_correlator");
      write(xml_seq_src, "backward_prop_corr", backward_prop_corr);
      pop(xml_seq_src);
      

      /*
       *  Write the sequential propagator component out to disk
       */
      {
	XMLBufferWriter file_xml;
	push(file_xml, "seqpropComponent");
	int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
	write(file_xml, "id", id);
	pop(file_xml);

	// Make a seqprop structure
	ChromaSeqProp_t seqprop_header;
	seqprop_header.invParam = input.param.invParam;
	seqprop_header.nonRelSeqProp  = input.param.nonRelSeqProp;
	seqprop_header.seq_src  = SeqSourceType(seq_src_value);
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
	  
	ostringstream outfile;
	outfile << input.prop.seqprop_files[seq_src_ctr] << "_component_s" << input.components[comp].spin << "_c" << input.components[comp].color;
 
	// Write the seqprop
	writeFermion(file_xml, record_xml, seq_prop_component,
		     outfile.str(), 
		     input.prop.seqprop_volfmt, QDPIO_SERIAL);
      }

	

    } /* End loop over components */
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

  END_CODE("seqpropComponent");

  exit(0);
}

