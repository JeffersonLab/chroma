// $Id: bar3ptfn.cc,v 1.22 2003-12-17 17:35:08 edwards Exp $
/*! \file
 * \brief Main program for computing 3pt functions
 *
 * Main program for computing 3pt functions
 */

#include "chroma.h"

using namespace QDP;

enum CfgType {
  CFG_TYPE_MILC,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_UNKNOWN
};

enum FermType {
  FERM_TYPE_WILSON,
  FERM_TYPE_UNKNOWN
};


/*
 * Input 
 */
struct IO_version_t
{
  int version;
};

// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType FermTypeP;

  multi1d<Real> Kappa;     // array of Wilson mass values

  CfgType cfg_type;        // storage order for stored gauge configuration
  int j_decay;             // direction to measure propagation

  bool Pt_src;             // point source
  bool Sl_src;             // shell source
  bool Pt_snk;             // point sink
  bool Sl_snk;             // shell sink

  int t_sink;

  multi1d<int> sink_mom;

  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.

  WvfKind Wvf_kind;        // Wave function kind: gauge invariant
  multi1d<Real> wvf_param; // Array of width's or other parameters
  //   for "shell" source/sink wave function
  multi1d<int> WvfIntPar;  // Array of iter numbers to approx. Gaussian or
  //   terminate CG inversion for Wuppertal smearing

  multi1d<int> Seq_src;

  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce;
};

struct Cfg_t
{
  string       cfg_file;
};

struct Bar3ptfn_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
};


// Reader for input parameters
void read(XMLReader& xml, const string& path, Bar3ptfn_input_t& input)
{
  XMLReader inputtop(xml, path);


  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  // Read in the IO_version
  try
  {
    read(inputtop, "IO_version/version", input.io_version.version);
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
    case 4 :
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

    {
      int input_Nc;
      read(paramtop, "Nc", input_Nc);
	
      if (input_Nc != Nc) {
	QDPIO::cerr << "Input parameter Nc=" << input_Nc \
		    <<  " different from qdp++ value." << endl;
	QDP_abort(1);
      }

      int input_Nd;
      read(paramtop, "Nd", input_Nd);

      if (input_Nd != Nd) {
	QDPIO::cerr << "Input parameter Nd=" << input_Nd \
		    << " different from qdp++ value." << endl;
	QDP_abort(1);
      }

      string ferm_type_str;
      read(paramtop, "FermTypeP", ferm_type_str);
      if (ferm_type_str == "WILSON") {
	input.param.FermTypeP = FERM_TYPE_WILSON;
      } else {
	input.param.FermTypeP = FERM_TYPE_UNKNOWN;
      }
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) 
    {
    case FERM_TYPE_WILSON :

      QDPIO::cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl;

      read(paramtop, "Kappa", input.param.Kappa);

      for (int i=0; i < input.param.Kappa.size(); ++i) {
	if (toBool(input.param.Kappa[i] < 0.0)) {
	  QDPIO::cerr << "Unreasonable value for Kappa." << endl;
	  QDPIO::cerr << "  Kappa[" << i << "] = " << input.param.Kappa[i] << endl;
	  QDP_abort(1);
	} else {
	  QDPIO::cout << " Spectroscopy Kappa: " << input.param.Kappa[i] << endl;
	}
      }

      break;

    default :
      QDPIO::cerr << "Fermion type not supported." << endl;
      if (input.param.FermTypeP == FERM_TYPE_UNKNOWN) {
	QDPIO::cerr << "  FermTypeP = UNKNOWN" << endl;
      }
      QDP_abort(1);
    }

    {
      string cfg_type_str;
      read(paramtop, "cfg_type", cfg_type_str);
      if (cfg_type_str == "SZIN") {
	input.param.cfg_type = CFG_TYPE_SZIN;
      } else {
	input.param.cfg_type = CFG_TYPE_UNKNOWN;
      }
    }

    read(paramtop, "j_decay", input.param.j_decay);
    if (input.param.j_decay < 0 || input.param.j_decay >= Nd) {
      QDPIO::cerr << "Bad value: j_decay = " << input.param.j_decay << endl;
      QDP_abort(1);
    }

    read(paramtop, "Pt_src", input.param.Pt_src);
    read(paramtop, "Sl_src", input.param.Sl_src);
    read(paramtop, "Pt_snk", input.param.Pt_snk);
    read(paramtop, "Sl_snk", input.param.Sl_snk);

    read(paramtop, "mom2_max", input.param.mom2_max);

    {
      string wvf_kind_str;
      read(paramtop, "Wvf_kind", wvf_kind_str);
      if (wvf_kind_str == "GAUGE_INV_GAUSSIAN") {
	input.param.Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN;
      } else {
	QDPIO::cerr << "Unsupported gauge-invariant Wvf_kind." << endl;
	QDPIO::cerr << "  Wvf_kind = " << wvf_kind_str << endl;
	QDP_abort(1);
      }
    }

    read(paramtop, "wvf_param", input.param.wvf_param);
    read(paramtop, "WvfIntPar", input.param.WvfIntPar);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);

    // Now we read in the information associated with the sequential sources
    read(paramtop, "t_sink", input.param.t_sink);
    read(paramtop, "Seq_src", input.param.Seq_src);
    read(paramtop, "sink_mom", input.param.sink_mom);

    for (int seq_src_ctr=0; seq_src_ctr<input.param.Seq_src.size(); ++seq_src_ctr) 
    {
      QDPIO::cout << "Computing sequential source of type "
		  << input.param.Seq_src[seq_src_ctr] << endl;
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg/cfg_file",input.cfg.cfg_file);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
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
  multi1d<int>      t_srce;
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
  multi1d<FormFac_Wilson_3Pt_fn_measurements_t>  bar;
};


// params
void write(BinaryWriter& bin, const Output_version_t& ver)
{
  write(bin, ver.out_version);
}

// params
void write(BinaryWriter& bin, const Param_t& param)
{
  write(bin, param.FermTypeP);
  write(bin, param.Kappa);

  write(bin, param.j_decay);
  write(bin, param.Pt_src);
  write(bin, param.Sl_src);
  write(bin, param.Pt_snk);
  write(bin, param.Sl_snk);
  write(bin, param.Wvf_kind);
  
  write(bin, param.t_sink);
  write(bin, param.sink_mom);

  write(bin, param.wvf_param);

  write(bin, param.WvfIntPar);
  write(bin, param.mom2_max);

  write(bin, param.Seq_src);
  write(bin, param.nrow);
}


// 
void write(BinaryWriter& bin, const FormFac_sequential_source_t& src)
{
  write(bin, src.seq_src_value);
  write(bin, src.t_srce);
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

  /*
   * Turn on the boundary conditions through the phase factors.
   *
   * NOTE: this is not an optimal solution: this factor stuff should be
   * set some other way
   */
  setph(input.param.boundary);              // initialize the BC factors


  // Sanity checks
  for (int i=0; i<Nd; ++i) {
    if (input.param.t_srce[i] < 0 || input.param.t_srce[i] >= input.param.nrow[i]) {
      QDPIO::cerr << "Quark propagator source coordinate incorrect." << endl;
      QDPIO::cerr << "t_srce[" << i << "] = " << input.param.t_srce[i] << endl;
      QDP_abort(1);
    }
  }

  if (input.param.t_sink < 0 || input.param.t_sink >= input.param.nrow[input.param.j_decay]) {
    QDPIO::cerr << "Sink time coordinate incorrect." << endl;
    QDPIO::cerr << "t_sink = " << input.param.t_sink << endl;
    QDP_abort(1);
  }

  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  // Check for unnecessary multiple occurances of kappas and/or wvf_params
  if (input.param.Kappa.size() > 1) {
    if (input.param.Sl_src == true) {
      for (int i=1; i < input.param.Kappa.size(); ++i) {
        for (int j=0; j<i; ++j) {
          if (toBool(input.param.Kappa[j] == input.param.Kappa[i])
              && toBool(input.param.wvf_param[j] == input.param.wvf_param[i])) {
            QDPIO::cerr << "Same kappa and wvf_param:" << endl;
            QDPIO::cerr << "  Kappa["     << i << "] = " << input.param.Kappa[i]     << endl;
            QDPIO::cerr << "  wvf_param[" << i << "] = " << input.param.wvf_param[i] << endl;
            QDPIO::cerr << "  Kappa["     << j << "] = " << input.param.Kappa[j]     << endl;
            QDPIO::cerr << "  wvf_param[" << j << "] = " << input.param.wvf_param[j] << endl;
            QDP_abort(1);
          }
        }
      }
    } else {
      for (int i=1; i < input.param.Kappa.size(); ++i) {
        for (int j=0; j<i; ++j) {
          if (toBool(input.param.Kappa[j] == input.param.Kappa[i])) {
            QDPIO::cerr  << "Same kappa without shell source or sink:" << endl;
            QDPIO::cerr << "  Kappa["     << i << "] = " << input.param.Kappa[i]     << endl;
            QDPIO::cerr << "  Kappa["     << j << "] = " << input.param.Kappa[j]     << endl;
            QDP_abort(1);
          }
        }
      }
    }
  }

  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  // Read in the configuration along with relevant information.
  QDPIO::cout << "Attempt to initialize the gauge field" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
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

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "IO_version");
  write(xml_out, "version", input.io_version.version);
  pop(xml_out);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 7);
  pop(xml_out);

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  xml_out.flush();

  // Big nested structure that is image of entire file
  Bar3ptfn_t  bar3pt;
  bar3pt.bar.resize(input.param.Kappa.size());
  bar3pt.output_version.out_version = 7;  // bump this up everytime something changes
  bar3pt.param = input.param; // copy entire structure


  XMLArrayWriter xml_array(xml_out, input.param.Kappa.size());
  push(xml_array, "Wilson_3Pt_fn_measurements");

  // Now loop over the various kappas
  for (int loop=0; loop < input.param.Kappa.size(); ++loop) 
  {
    QDPIO::cout << "Mass loop = " << loop << endl;
  
    push(xml_array);
    Write(xml_array, loop);

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
   

    // Big nested structure that is image of all form-factors
//    FormFac_Wilson_3Pt_fn_measurements_t  formfacs;
    bar3pt.bar[loop].output_version = 1;  // bump this up everytime something changes
    bar3pt.bar[loop].seqsrc.resize(input.param.Seq_src.size());

    XMLArrayWriter  xml_seq_src(xml_array, input.param.Seq_src.size());
    push(xml_seq_src, "Sequential_source");

    for (int seq_src_ctr = 0; seq_src_ctr < input.param.Seq_src.size(); ++seq_src_ctr) 
    {
      push(xml_seq_src);
      Write(xml_seq_src, seq_src_ctr);

      LatticePropagator seq_quark_prop;
      int seq_src_value = input.param.Seq_src[seq_src_ctr];

      // Read the sequential propagator
      QDPIO::cout << "Attempt to read backward propagator" << endl;
      XMLReader seqprop_xml;
      {
        stringstream prop_file;
        prop_file << "seqprop_" << loop << "_" << seq_src_value;
        readSzinQprop(seqprop_xml, seq_quark_prop, prop_file.str());

	write(xml_seq_src, "Backward_prop_info", seqprop_xml);
      }

      QDPIO::cout << "Backward propagator successfully read" << endl;
   
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
      write(xml_seq_src, "t_srce", input.param.t_srce);
      write(xml_seq_src, "t_sink", input.param.t_sink);
      write(xml_seq_src, "sink_mom", input.param.sink_mom);
	
      bar3pt.bar[loop].seqsrc[seq_src_ctr].seq_src_value = seq_src_value;
      bar3pt.bar[loop].seqsrc[seq_src_ctr].t_srce        = input.param.t_srce;
      bar3pt.bar[loop].seqsrc[seq_src_ctr].t_sink        = input.param.t_sink;
      bar3pt.bar[loop].seqsrc[seq_src_ctr].sink_mom      = input.param.sink_mom;
	
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
      LatticePropagator seq_quark_prop_tmp = seq_quark_prop;

      if (input.param.Sl_src) {
        sink_smear2(u, seq_quark_prop_tmp, input.param.Wvf_kind, input.param.wvf_param[loop],
                    input.param.WvfIntPar[loop], input.param.j_decay);
      }

      // Compute the 2pt function at the source - used in later normalizations
      Complex seq_hadron_0 =
        peekSite(trace(adj(Gamma(G5)*seq_quark_prop_tmp*Gamma(G5))), input.param.t_srce);

      Write(xml_seq_src, seq_hadron_0);
      bar3pt.bar[loop].seqsrc[seq_src_ctr].seq_hadron_0 = seq_hadron_0;

      // Now the 3pt contractions
      SftMom phases(input.param.mom2_max, input.param.sink_mom, false, input.param.j_decay);
      FormFac(bar3pt.bar[loop].seqsrc[seq_src_ctr].formFacs, 
	      u, quark_propagator, seq_quark_prop, phases, input.param.t_srce[input.param.j_decay]);

      pop(xml_seq_src);   // elem
    } // end loop over sequential sources

    pop(xml_seq_src);  // Sequential_source

//    BinaryWriter  bin_out("bar3ptfn.dat");
//    write(bin_out, bar3ptfn);
//    bin_out.close();

    pop(xml_array);     // elem
  } // end loop over the kappa value

  pop(xml_array);  // Wilson_3Pt_fn_measurements

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
