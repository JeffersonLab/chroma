// $Id: spectrum_w.cc,v 1.11 2003-10-01 20:23:46 edwards Exp $
//
//! \file
//  \brief Main code for propagator generation
//
//  $Log: spectrum_w.cc,v $
//  Revision 1.11  2003-10-01 20:23:46  edwards
//  Now supports baryons and currents. Changed reading to use a
//  structure. Pushed all control vars to this structure.
//
//  Revision 1.10  2003/09/10 18:04:22  edwards
//  Changed to new form of XMLReader - a clone.
//
//  Revision 1.9  2003/09/02 15:52:02  edwards
//  Added return 0 at end of main. Some kind of return is required in C++.
//
//  Revision 1.8  2003/08/27 22:08:41  edwards
//  Start major push to using xml.
//
//  Revision 1.7  2003/08/27 20:05:20  edwards
//  Removed use of seed. Do not need random numbers here.
//
//  Revision 1.6  2003/06/24 03:25:06  edwards
//  Changed from nml to xml.
//
//  Revision 1.5  2003/06/08 05:00:25  edwards
//  Added some flush to nml_out.
//
//  Revision 1.4  2003/05/22 17:35:36  flemingg
//  Added stripper for spectrum_w.  Also, minor change to spectrum_w.cc
//  to make it compatible with the stripper.
//
//  Revision 1.3  2003/05/13 22:00:50  flemingg
//  I'm done with spectrum_w and the test files for now. I'm happy
//  enough with the output format. Somebody please write a stripper.
//
//  Revision 1.2  2003/05/08 23:02:12  flemingg
//  Initial version of spectrum_w.  It compiles and reproduces szin meson
//  correlation functions up to a sign.  Namelist input and output starting
//  to evolve away from szin model: a work in progress.  Code runs but
//  more testing needed including resolving the cause of the different signs.
//

#include <iostream>
#include <cstdio>

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
  FermType     FermTypeP;
  int          numKappa;   // number of Wilson masses
  multi1d<Real>  Kappa;    // array of Wilson mass values - someday use a mass instead

  CfgType cfg_type;        // storage order for stored gauge configuration
  int j_decay;             // direction to measure propagation

  bool Pt_src;             // point source
  bool Sl_src;             // shell source
  bool Pt_snk;             // point sink
  bool Sl_snk;             // shell sink

  bool MesonP;             // Meson spectroscopy
  bool CurrentP;           // Meson currents
  bool BaryonP;            // Baryons spectroscopy

  bool time_rev;           // Use time reversal in baryon spectroscopy


  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  bool avg_equiv_mom;      // average over equivalent momenta
  WvfKind Wvf_kind;        // Wave function kind: gauge invariant
  multi1d<Real> wvf_param; // Array of width's or other parameters
  //   for "shell" source/sink wave function
  multi1d<int> WvfIntPar;  // Array of iter numbers to approx. Gaussian or
  //   terminate CG inversion for Wuppertal smearing

  multi1d<int> disk_prop;
  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce;
};

struct Cfg_t
{
  string       cfg_file;
};

struct Input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
};


// Reader for input parameters
void read(XMLReader& xml, const string& path, Input_t& input)
{
  XMLReader inputtop(xml, path);

  // Defaults
  input.param.time_rev = false;
  input.param.CurrentP = false;
  input.param.BaryonP  = false;


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
    cerr << "Error reading data: " << e << endl;
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
    case 6 :
      /**************************************************************************/
      break;

    case 7:
      /**************************************************************************/
      read(paramtop, "CurrentP", input.param.CurrentP);
      read(paramtop, "BaryonP", input.param.BaryonP);
      break;

    default :
      /**************************************************************************/

      cerr << "Input parameter version " << input.io_version.version << " unsupported." << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    cerr << "Error reading data: " << e << endl;
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
	cerr << "Input parameter Nc=" << input_Nc \
	     <<  " different from qdp++ value." << endl;
	QDP_abort(1);
      }

      int input_Nd;
      read(paramtop, "Nd", input_Nd);

      if (input_Nd != Nd) {
	cerr << "Input parameter Nd=" << input_Nd \
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
    switch (input.param.FermTypeP) {
    case FERM_TYPE_WILSON :

//	cout << " SPECTRUM_W: Spectroscopy for Wilson fermions" << endl;
      QDP_info(" SPECTRUM_W: Spectroscopy for Wilson fermions");

      read(paramtop, "numKappa", input.param.numKappa);
      read(paramtop, "Kappa", input.param.Kappa);

      for (int i=0; i < input.param.numKappa; ++i) {
	if (toBool(input.param.Kappa[i] < 0.0)) {
	  cerr << "Unreasonable value for Kappa." << endl;
	  cerr << "  Kappa[" << i << "] = " << input.param.Kappa[i] << endl;
	  QDP_abort(1);
	} else {
	  cout << " Spectroscopy Kappa: " << input.param.Kappa[i] << endl;
	}
      }

      break;

    default :
      cerr << "Fermion type not supported." << endl;
      if (input.param.FermTypeP == FERM_TYPE_UNKNOWN) {
	cerr << "  FermTypeP = UNKNOWN" << endl;
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
      cerr << "Bad value: j_decay = " << input.param.j_decay << endl;
      QDP_abort(1);
    }

    read(paramtop, "Pt_src", input.param.Pt_src);
    read(paramtop, "Sl_src", input.param.Sl_src);
    read(paramtop, "Pt_snk", input.param.Pt_snk);
    read(paramtop, "Sl_snk", input.param.Sl_snk);

    read(paramtop, "MesonP", input.param.MesonP);

    read(paramtop, "mom2_max", input.param.mom2_max);
    read(paramtop, "avg_equiv_mom", input.param.avg_equiv_mom);

    {
      string wvf_kind_str;
      read(paramtop, "Wvf_kind", wvf_kind_str);
      if (wvf_kind_str == "GAUGE_INV_GAUSSIAN") {
	input.param.Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN;
      } else {
	cerr << "Unsupported gauge-invariant Wvf_kind." << endl;
	cerr << "  Wvf_kind = " << wvf_kind_str << endl;
	QDP_abort(1);
      }
    }

    read(paramtop, "wvf_param", input.param.wvf_param);
    read(paramtop, "WvfIntPar", input.param.WvfIntPar);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);
  }
  catch (const string& e) 
  {
    cerr << "Error reading data: " << e << endl;
    throw;
  }



  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg/cfg_file",input.cfg.cfg_file);
  }
  catch (const string& e) 
  {
    cerr << "Error reading data: " << e << endl;
    throw;
  }
}


//
// Main program
//
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  Input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/Input", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  /*
   * Sanity checks
   */
  // Check that only one type of source smearing is specified
  // Must match how the stored propagator was smeared
  if (input.param.Pt_src && input.param.Sl_src) {
    cerr << "Error if Pt_src and Sl_src are both set." << endl;
    cerr << "Choose the one which matches the stored propagators." << endl;
    QDP_abort(1);
  }

  // Figure out what to do about boundary conditions
  // GTF HACK: only allow periodic boundary conditions
  for (int i=0; i<Nd; ++i) {
    if (input.param.boundary[i] != 1) {
      cerr << "Only periodic boundary conditions supported." << endl;
      cerr << "  boundary[" << i << "] = " << input.param.boundary[i] << endl;
      QDP_abort(1);
    }
  }

  for (int i=0; i<Nd; ++i) {
    if (input.param.t_srce[i] < 0 || input.param.t_srce[i] >= input.param.nrow[i]) {
      cerr << "Quark propagator source coordinate incorrect." << endl;
      cerr << "t_srce[" << i << "] = " << input.param.t_srce[i] << endl;
      QDP_abort(1);
    }
  }


  cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  multi1d<LatticeColorMatrix> u(Nd);

  cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    cout << " x " << input.param.nrow[i];
  }
  cout << endl;

  // Read in the configuration along with relevant information.
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
  push(xml_out, "spectrum_w");

  xml_out << xml_in;     // save a copy of the input

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 6);
  pop(xml_out);

  xml_out.flush();


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

  // Next check the gauge field configuration by reunitarizing.
  {
    LatticeBoolean lbad;
    int numbad;
    multi1d<LatticeColorMatrix> u_tmp(Nd);
    u_tmp = u;
    int mu;
    for (mu=0; mu < Nd; ++mu) {
      reunit(u_tmp[mu], lbad, numbad, REUNITARIZE_ERROR);
    }
  }

  // Initialize the slow Fourier transform phases
  SftMom phases(input.param.mom2_max, input.param.avg_equiv_mom, input.param.j_decay);

  // Keep an array of all the xml output buffers
  XMLArrayWriter xml_array(xml_out,input.param.numKappa);
  push(xml_array, "Wilson_hadron_measurements");


  // Flags
  int t0      = input.param.t_srce[input.param.j_decay];
  int bc_spec = input.param.boundary[input.param.j_decay];

  // Now loop over the various fermion masses
  for (int loop=0; loop < input.param.numKappa; ++loop)
  {
    // Read the quark propagator
    LatticePropagator quark_propagator;
    {
      XMLReader prop_xml;
      stringstream prop_file;
      prop_file << "propagator_" << loop;
      readSzinQprop(prop_xml, quark_propagator, prop_file.str());
    }

    push(xml_array);         // next array element - name auto-written
    Write(xml_array, loop);
    write(xml_array, "Kappa_mes", input.param.Kappa[loop]);
    write(xml_array, "t_srce", input.param.t_srce);

    // Do the mesons first
    if (input.param.MesonP) 
    {
      // Construct {Point|Shell}-Point mesons, if desired
      if (input.param.Pt_snk) 
      {
	if (input.param.Pt_src)
	  mesons(quark_propagator, quark_propagator, phases, t0,
		 xml_array, "Point_Point_Wilson_Mesons");
        
	if (input.param.Sl_src)
	  mesons(quark_propagator, quark_propagator, phases, t0,
		 xml_array, "Shell_Point_Wilson_Mesons");
      } // end if (Pt_snk)

      // Convolute the quark propagator with the sink smearing function.
      // Make a copy of the quark propagator and then overwrite it with
      // the convolution. 
      if (input.param.Sl_snk) 
      {
	LatticePropagator quark_prop_smr;
	quark_prop_smr = quark_propagator;
	sink_smear2(u, quark_prop_smr, input.param.Wvf_kind, input.param.wvf_param[loop],
		    input.param.WvfIntPar[loop], input.param.j_decay);

	if (input.param.Pt_src)
	  mesons(quark_prop_smr, quark_prop_smr, phases, t0,
		 xml_array, "Point_Shell_Wilson_Mesons");

	if (input.param.Sl_src)
	  mesons(quark_prop_smr, quark_prop_smr, phases, t0,
		 xml_array, "Shell_Shell_Wilson_Mesons");
      } // end if (Sl_snk)

    } // end if (MesonP)


    // Do the currents next
    if (input.param.CurrentP) 
    {
      // Construct the rho vector-current and the pion axial current divergence
      if (input.param.Pt_src)
	curcor2(u, quark_propagator, quark_propagator, phases, 
		t0, input.param.j_decay, 4,
		xml_array, "Point_Point_Meson_Currents");
        
      if (input.param.Sl_src)
	curcor2(u, quark_propagator, quark_propagator, phases, 
		t0, input.param.j_decay, 4,
		xml_array, "Shell_Point_Meson_Currents");
    } // end if (CurrentP)


    // Do the baryons
    if (input.param.BaryonP) 
    {
      // Construct {Point|Shell}-Point mesons, if desired
      if (input.param.Pt_snk) 
      {
	if (input.param.Pt_src)
	  baryon(quark_propagator, phases, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Point_Point_Wilson_Baryons");
        
	if (input.param.Sl_src)
	  baryon(quark_propagator, phases, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Shell_Point_Wilson_Baryons");
      } // end if (Pt_snk)

      // Convolute the quark propagator with the sink smearing function.
      // Make a copy of the quark propagator and then overwrite it with
      // the convolution. 
      if (input.param.Sl_snk) 
      {
	LatticePropagator quark_prop_smr;
	quark_prop_smr = quark_propagator;
	sink_smear2(u, quark_prop_smr, input.param.Wvf_kind, input.param.wvf_param[loop],
		    input.param.WvfIntPar[loop], input.param.j_decay);
	if (input.param.Pt_src)
	  baryon(quark_propagator, phases, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Point_Shell_Wilson_Baryons");
        
	if (input.param.Sl_src)
	  baryon(quark_propagator, phases, 
		 t0, bc_spec, input.param.time_rev, 	
		 xml_array, "Shell_Shell_Wilson_Baryons");
      } // end if (Sl_snk)

    } // end if (BaryonP)

    pop(xml_array);  // array element

  } // end for(loop)

  pop(xml_array);  // Wilson_spectroscopy
  pop(xml_out);  // spectrum_w

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
