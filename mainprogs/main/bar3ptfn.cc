// $Id: bar3ptfn.cc,v 1.20 2003-10-09 20:32:37 edwards Exp $
//
// $Log: bar3ptfn.cc,v $
// Revision 1.20  2003-10-09 20:32:37  edwards
// Changed all cout/cerr to QDPIO::cout/cerr. Change QDP_info calls
// to use QDPIO::cout.
//
// Revision 1.19  2003/10/06 21:15:46  edwards
// Use new input format.
//
// Revision 1.18  2003/09/11 15:27:18  edwards
// Turned on some flush's.
//
// Revision 1.17  2003/09/11 15:25:32  edwards
// Added some diagnostic output.
//
// Revision 1.16  2003/09/11 15:17:34  edwards
// Added try/catch to DATA readers.
//
// Revision 1.15  2003/09/10 18:04:22  edwards
// Changed to new form of XMLReader - a clone.
//
// Revision 1.14  2003/08/27 22:08:41  edwards
// Start major push to using xml.
//
// Revision 1.13  2003/08/27 20:04:14  edwards
// Changed readSzin to return an xml header.
//
// Revision 1.12  2003/07/04 17:08:36  edwards
// Added more Seq_src types.
//
// Revision 1.11  2003/06/25 16:12:04  edwards
// Changed from nml to xml.
//
// Revision 1.10  2003/06/08 05:02:27  edwards
// Added some flush to nml_out.
//
// Revision 1.9  2003/05/30 02:37:40  flemingg
// A message printed to stdout was printing the wrong thing cut and pasted
// from spectrum_w.cc.  Replaced with the intended thing from szin bar3ptfn.m
//
// Revision 1.8  2003/05/14 06:07:03  flemingg
// Should be done playing around with the input and output formats
// for bar3ptfn with this version.
//
//

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

  int numKappa;            // number of Wilson masses
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

  int numSeq_src;
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
      QDPIO::cerr << "Input parameter version " << input.io_version.version << " unsupported." << endl;
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

      read(paramtop, "numKappa", input.param.numKappa);
      read(paramtop, "Kappa", input.param.Kappa);

      for (int i=0; i < input.param.numKappa; ++i) {
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

    read(paramtop, "t_sink", input.param.t_sink);

    // Now we read in the information associated with the sequential sources
    read(paramtop, "numSeq_src", input.param.numSeq_src);

    read(paramtop, "sink_mom", input.param.sink_mom);

    // Now read in the particular Sequential Sources we are evaluating
    read(paramtop, "Seq_src", input.param.Seq_src);

    for (int seq_src_ctr=0; seq_src_ctr<input.param.numSeq_src; ++seq_src_ctr) 
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

  // Figure out what to do about boundary conditions
  // GTF HACK: only allow periodic boundary conditions
  for (int i=0; i<Nd; ++i) {
    if (input.param.boundary[i] != 1) {
      QDPIO::cerr << "Only periodic input.param.boundary conditions supported." << endl;
      QDPIO::cerr << "  input.param.boundary[" << i << "] = " << input.param.boundary[i] << endl;
      QDP_abort(1);
    }
  }

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
  if (input.param.numKappa > 1) {
    if (input.param.Sl_src == true) {
      for (int i=1; i < input.param.numKappa; ++i) {
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
      for (int i=1; i < input.param.numKappa; ++i) {
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
  write(xml_out, "out_version", 6);
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

  XMLArrayWriter xml_array(xml_out, input.param.numKappa);
  push(xml_array, "Wilson_3Pt_fn_measurements");

  // Now loop over the various kappas
  for (int loop=0; loop < input.param.numKappa; ++loop) 
  {
    QDPIO::cout << "Mass loop = " << loop << endl;
  
    push(xml_array);
    Write(xml_array, loop);

    // Read the quark propagator
    QDPIO::cout << "Attempt to read forward propagator" << endl;
  
    LatticePropagator quark_propagator;
    {
      XMLReader prop_xml;
      stringstream prop_file;
      prop_file << "propagator_" << loop;
      readSzinQprop(prop_xml, quark_propagator, prop_file.str());

      write(xml_array, "Forward_prop_info", prop_xml);
    }

    QDPIO::cout << "Forward propagator successfully read" << endl;
   

    XMLArrayWriter  xml_seq_src(xml_array, input.param.numSeq_src);
    push(xml_seq_src, "Sequential_source");

    for (int seq_src_ctr = 0; seq_src_ctr < input.param.numSeq_src; ++seq_src_ctr) 
    {
      push(xml_seq_src);
      Write(xml_seq_src, seq_src_ctr);

      LatticePropagator seq_quark_prop;
      int seq_src_value = input.param.Seq_src[seq_src_ctr];

      // Read the sequential propagator
      QDPIO::cout << "Attempt to read backward propagator" << endl;
   
      {
	XMLReader seqprop_xml;
        stringstream prop_file;
        prop_file << "seqprop_" << loop << "_" << seq_src_value;
        readSzinQprop(seqprop_xml, seq_quark_prop, prop_file.str());

	write(xml_array, "Backward_prop_info", seqprop_xml);
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

//      xml_seq_src.flush();

      // Now the 3pt contractions
      SftMom phases(input.param.mom2_max, input.param.sink_mom, false, input.param.j_decay);
      FormFac(u, quark_propagator, seq_quark_prop, phases, input.param.t_srce[input.param.j_decay],
              xml_seq_src);

      pop(xml_seq_src);   // elem
//      xml_seq_src.flush();
    } // end loop over sequential sources

    pop(xml_seq_src);  // Sequential_source
    pop(xml_array);     // elem
  } // end loop over the kappa value

  pop(xml_array);  // Wilson_3Pt_fn_measurements

  // Close the namelist output file NMLDAT
  pop(xml_out);     // bar3ptfn

  xml_in.close();
  xml_out.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
