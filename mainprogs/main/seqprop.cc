// $Id: seqprop.cc,v 1.8 2004-01-07 04:44:26 edwards Exp $
/*! \file
 *  \brief Main code for sequential propagator generation
 */

// #error "CODE NOT READY YET"


#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;


/*
 * Input 
 */

//! Parameters for sources and sinks
struct SrceSinkParam_t
{
  bool             Pt_src;   // point source
  bool             Sl_src;   // shell source
  bool             Pt_snk;   // point sink
  bool             Sl_snk;   // shell sink

  SmearingParam_t  smearParam;

  multi1d<int>     t_srce;
  multi1d<int>     sink_mom;
  int              t_sink;
};

// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType         FermTypeP;
  FermActType      FermAct;
  multi1d<Real>    Mass;       // Quark Mass and **NOT** kappa
 
  ChiralParam_t    chiralParam;
  SrceSinkParam_t  srceSinkParam;
  InvertParam_t    invParam;

  CfgType          cfg_type;   // storage order for stored gauge configuration
  PropType         prop_type;  // storage order for stored propagator

  int              j_decay;    // decay direction

  multi1d<int>     seq_src;    // integer array holding sequential source numbers

  multi1d<int>     nrow;
  multi1d<int>     boundary;
};

struct Seqprop_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
};


//
void read(XMLReader& xml, const string& path, SrceSinkParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "Pt_src", param.Pt_src);
  read(paramtop, "Sl_src", param.Sl_src);
  read(paramtop, "Pt_snk", param.Pt_snk);
  read(paramtop, "Sl_snk", param.Sl_snk);

  if (param.Sl_src || param.Sl_snk)
    read(paramtop, "smearParam", param.smearParam);

  read(paramtop, "t_srce", param.t_srce);
  read(paramtop, "t_sink", param.t_sink);
  read(paramtop, "sink_mom", param.sink_mom);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, Seqprop_input_t& input)
{
  XMLReader inputtop(xml, path);

  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  // Read in the IO_version
  try
  {
    read(inputtop, "IO_version", input.io_version);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Currently, in the supported IO versions, there is only a small difference
  // in the inputs. So, to make code simpler, extract the common bits 

  // Read the uncommon bits first
  try
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    switch (input.io_version.version) 
    {
      /**************************************************************************/
    case 1 :
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << input.io_version.version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read the common bits
  try 
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    read(paramtop, "FermTypeP", input.param.FermTypeP);
    read(paramtop, "FermAct", input.param.FermAct);

    if (paramtop.count("Mass") != 0)
    {
      read(paramtop, "Mass", input.param.Mass);

      if (paramtop.count("Kappa") != 0)
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      
      cout << "conver kappa to mass" << endl;
      input.param.Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }    

#if 0
    for (int i=0; i < input.param.Mass.size(); ++i) {
      if (toBool(input.param.Mass[i] < 0.0)) {
	QDPIO::cerr << "Unreasonable value for Mass." << endl;
	QDPIO::cerr << "  Mass[" << i << "] = " << input.param.Mass[i] << endl;
	QDP_abort(1);
      } else {
	QDPIO::cout << " Spectroscopy Mass: " << input.param.Mass[i] << endl;
      }
    }
#endif

    if (paramtop.count("ChiralParam") != 0)
      read(paramtop, "ChiralParam", input.param.chiralParam);

    read(paramtop, "cfg_type", input.param.cfg_type);
    read(paramtop, "prop_type", input.param.prop_type);

    read(paramtop, "SrceSinkParam", input.param.srceSinkParam);
    read(paramtop, "InvertParam", input.param.invParam);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg", input.cfg);
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
  for (int i=0; i<Nd; ++i) {
    if (input.param.srceSinkParam.t_srce[i] < 0 || input.param.srceSinkParam.t_srce[i] >= input.param.nrow[i]) {
      QDPIO::cerr << "Quark propagator source coordinate incorrect." << endl;
      QDPIO::cerr << "t_srce[" << i << "] = " << input.param.srceSinkParam.t_srce[i] << endl;
      QDP_abort(1);
    }
  }

  if (input.param.srceSinkParam.t_sink < 0 || input.param.srceSinkParam.t_sink >= input.param.nrow[input.param.j_decay]) {
    QDPIO::cerr << "Sink time coordinate incorrect." << endl;
    QDPIO::cerr << "t_sink = " << input.param.srceSinkParam.t_sink << endl;
    QDP_abort(1);
  }

  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  // Check for unnecessary multiple occurances of masss and/or wvf_params
  if (input.param.Mass.size() > 1) {
    if (input.param.srceSinkParam.Sl_src == true) {
      for (int i=1; i < input.param.Mass.size(); ++i) {
        for (int j=0; j<i; ++j) {
          if (toBool(input.param.Mass[j] == input.param.Mass[i])
              && toBool(input.param.srceSinkParam.smearParam.wvf_param[j] == input.param.srceSinkParam.smearParam.wvf_param[i])) {
            QDPIO::cerr << "Same Mass and wvf_param:" << endl;
            QDPIO::cerr << "  Mass["     << i << "] = " << input.param.Mass[i]     << endl;
            QDPIO::cerr << "  wvf_param[" << i << "] = " << input.param.srceSinkParam.smearParam.wvf_param[i] << endl;
            QDPIO::cerr << "  Mass["     << j << "] = " << input.param.Mass[j]     << endl;
            QDPIO::cerr << "  wvf_param[" << j << "] = " << input.param.srceSinkParam.smearParam.wvf_param[j] << endl;
            QDP_abort(1);
          }
        }
      }
    } else {
      for (int i=1; i < input.param.Mass.size(); ++i) {
        for (int j=0; j<i; ++j) {
          if (toBool(input.param.Mass[j] == input.param.Mass[i])) {
            QDPIO::cerr  << "Same Mass without shell source or sink:" << endl;
            QDPIO::cerr << "  Mass["     << i << "] = " << input.param.Mass[i]     << endl;
            QDPIO::cerr << "  Mass["     << j << "] = " << input.param.Mass[j]     << endl;
            QDP_abort(1);
          }
        }
      }
    }
  }

  QDPIO::cout << "\n     Gauge group: SU(" << Nc << ")" << endl;

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
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "seqprop");

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
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  xml_out.flush();


  //------------------ Start main body of calculations -----------------------------
  int ncg_had = 0;			/* Initialise iteration counter */

  /*
   * Construct fermionic BC. Need one for LatticeFermion and multi1d<LatticeFermion>
   * Note, the handle is on an ABSTRACT type
   */
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(input.param.boundary));

  /*
   *  Now loop over the various masses
   */
  XMLArrayWriter xml_array(xml_out, input.param.Mass.size());
  push(xml_array, "Seqprop_calc");
    
  for(int loop = 0; loop < input.param.Mass.size(); loop++)
  {
    QDPIO::cout << "Mass loop = " << loop << endl;
  
    push(xml_array);
    Write(xml_array, loop);

    //
    // Initialize fermion action
    //
    Real Mass_meas = input.param.Mass[loop];

#if 1
    UnprecWilsonFermAct S_f(fbc,Mass_meas);
#else
    UnprecDWFermActArray S_f(fbc_a,
			     input.param.chiralParam.OverMass, 
			     Mass_meas, 
			     input.param.chiralParam.N5);
//  UnprecDWFermAct S_f(fbc_a, WilsonMass, m_q);
#endif

    // Read the quark propagator
    QDPIO::cout << "Attempt to read forward propagator" << endl;
  
    LatticePropagator quark_propagator;
    XMLReader prop_xml;
    {
      stringstream prop_file;
      prop_file << "propagator_" << loop;
      readSzinQprop(prop_xml, quark_propagator, prop_file.str());

      write(xml_array, "Forward_prop_info", prop_xml);
    }

    QDPIO::cout << "Forward propagator successfully read" << endl;
   

    if (input.param.srceSinkParam.Sl_snk)
    {
      // Do the sink smearing BEFORE the interpolating operator
      sink_smear2(u, quark_propagator, 
		  input.param.srceSinkParam.smearParam.Wvf_kind, 
		  input.param.srceSinkParam.smearParam.wvf_param[loop], 
		  input.param.srceSinkParam.smearParam.WvfIntPar[loop], 
		  input.param.j_decay);
    }

    for(int seq_src_ctr = 0; seq_src_ctr < input.param.seq_src.size(); seq_src_ctr++)
    {
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
		     input.param.srceSinkParam.t_sink, 
		     input.param.srceSinkParam.sink_mom, 
		     input.param.j_decay, 
		     seq_src_value);
      }
      else if ((10 <= seq_src_value) && (seq_src_value <= 20))
      {
	// Computation of the Meson sequential source
	mesonSeqSource(quark_propagator, quark_prop_src, 
		       input.param.srceSinkParam.t_sink, 
		       input.param.srceSinkParam.sink_mom, 
		       input.param.j_decay, 
		       seq_src_value);
      }
      else{
	QDP_error_exit("Unknown sequential source type", seq_src_value);
      }

      if (input.param.srceSinkParam.Sl_snk)
      {
	// Do the sink smearing AFTER the interpolating operator
	sink_smear2(u, quark_prop_src, 
		    input.param.srceSinkParam.smearParam.Wvf_kind, 
		    input.param.srceSinkParam.smearParam.wvf_param[loop], 
		    input.param.srceSinkParam.smearParam.WvfIntPar[loop], 
		    input.param.j_decay);
      }

      /*
       *  Compute the full propagator.
       */
      LatticePropagator seq_quark_prop = zero;
      XMLBufferWriter xml_buf;
      int ncg_had;
      {
	Handle<const ConnectState> state(S_f.createState(u));  // inserts any BC

	quarkProp4(seq_quark_prop, xml_buf, quark_prop_src,
		   S_f, state, 
		   input.param.invParam.invType, 
		   input.param.invParam.RsdCG, 
		   input.param.invParam.MaxCG, 
		   ncg_had);
      }

      xml_out << xml_buf;

      /*
       *  Write the sequential propagator out to disk
       *
       *  We need to dump some sort of header.  At the very least, we
       *  should dump the Mass values
       */

      
      /*
       *  Now create the name of the propagator value
       */
      {
	stringstream seqprop_file;
	seqprop_file << "seqprop_" << loop << "_seq_src_value";
	writeSzinQprop(seq_quark_prop, seqprop_file.str(), Mass_meas);
      }

      /*
       *  In the case of the pion, we know that the exponentiated propagator
       *  back to the source should be the pion correlator at time-slice
       *  zero, and so will write this out
       */

      if(seq_src_value == 10)
      {
	Complex pion_src;
	seqPionTest(pion_src, seq_quark_prop, input.param.srceSinkParam.t_srce);
	
	push(xml_out,"Seq_propagator_test");
	Write(xml_out, pion_src);
	pop(xml_out);
      }
    } /* end loop over sequential sources */
      
  } /* end loop over the mass value */

  pop(xml_array);  // Seqprop_calc

  push(xml_out,"Relaxation_Iterations");
  Write(xml_out, ncg_had);
  pop(xml_out);

  pop(xml_out);

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE("seqprop");

  exit(0);
}

