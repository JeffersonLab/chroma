// $Id: spectrum_w.cc,v 1.8 2003-08-27 22:08:41 edwards Exp $
//
//! \file
//  \brief Main code for propagator generation
//
//  $Log: spectrum_w.cc,v $
//  Revision 1.8  2003-08-27 22:08:41  edwards
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
} ;

enum FermType {
  FERM_TYPE_WILSON,
  FERM_TYPE_UNKNOWN
} ;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  int i, j ;

  // Parameters which must be determined from the namelist input
  // and written to the namelist output

  int version ;              // input parameter version

  FermType FermTypeP ;
  // GTF GRIPE: I would prefer masses rather than kappa values here,
  //   but I'll save that for another day.
  int numKappa ;            // number of Wilson masses
  multi1d<Real> Kappa ;     // array of Wilson mass values

  CfgType cfg_type ;        // storage order for stored gauge configuration
  int j_decay ;             // direction to measure propagation

  bool Pt_src ;             // point source
  bool Sl_src ;             // shell source
  bool Pt_snk ;             // point sink
  bool Sl_snk ;             // shell sink

  bool MesonP ;             // Meson spectroscopy

  int mom2_max ;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  bool avg_equiv_mom ;      // average over equivalent momenta
  WvfKind Wvf_kind ;        // Wave function kind: gauge invariant
  multi1d<Real> wvf_param ; // Array of width's or other parameters
                            //   for "shell" source/sink wave function
  multi1d<int> WvfIntPar ;  // Array of iter numbers to approx. Gaussian or
                            //   terminate CG inversion for Wuppertal smearing

  multi1d<int> disk_prop ;
  multi1d<int> nrow(Nd) ;
  multi1d<int> boundary(Nd) ;
  multi1d<int> t_srce(Nd) ;

  string cfg_file ;

  // Instantiate namelist reader for DATA
  XMLReader xml_in("DATA") ;
  string xml_in_root = "/spectrum_w";

  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  string path = xml_in_root + "/IO_version";
  ReadPath(xml_in, path, version) ;

  switch (version) {

  /**************************************************************************/
  case 5 :
  /**************************************************************************/

    path = xml_in_root + "/param"; // push into 'param' group

    {
      int input_Nc ;
      read(xml_in, path + "/Nc", input_Nc) ;

      if (input_Nc != Nc) {
        cerr << "Input parameter Nc=" << input_Nc \
          <<  " different from qdp++ value." << endl ;
        QDP_abort(1) ;
      }

      int input_Nd ;
      read(xml_in, path + "/Nd", input_Nd) ;

      if (input_Nd != Nd) {
        cerr << "Input parameter Nd=" << input_Nd \
          << " different from qdp++ value." << endl ;
        QDP_abort(1) ;
      }

      int ferm_type_int ;
      read(xml_in, path + "/FermTypeP", ferm_type_int) ;
      switch (ferm_type_int) {
      case 1 :
        FermTypeP = FERM_TYPE_WILSON ;
        break ;
      default :
        FermTypeP = FERM_TYPE_UNKNOWN ;
      }
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (FermTypeP) {
    case FERM_TYPE_WILSON :

      cout << " SPECTRUM_W: Spectroscopy for Wilson fermions" << endl ;

      ReadPath(xml_in, path, numKappa) ;

      Kappa.resize(numKappa) ;
      ReadPath(xml_in, path, Kappa) ;

      for (i=0; i < numKappa; ++i) {
        if (toBool(Kappa[i] < 0.0)) {
          cerr << "Unreasonable value for Kappa." << endl ;
          cerr << "  Kappa[" << i << "] = " << Kappa[i] << endl ;
          QDP_abort(1) ;
        } else {
          cout << " Spectroscopy Kappa: " << Kappa[i] << endl ;
        }
      }

      break ;

    default :
      cerr << "Fermion type not supported." << endl ;
      if (FermTypeP == FERM_TYPE_UNKNOWN) {
        cerr << "  FermTypeP = UNKNOWN" << endl ;
      }
      QDP_abort(1) ;
    }

    {
      int input_cfg_type ;
      read(xml_in, path + "/cfg_type", input_cfg_type) ;
      switch (input_cfg_type) {
      case 1 :
        cfg_type = CFG_TYPE_SZIN ;
        break ;
      default :
        cfg_type = CFG_TYPE_UNKNOWN ;
      }
    }

    ReadPath(xml_in, path, j_decay) ;
    if (j_decay < 0 || j_decay >= Nd) {
      cerr << "Bad value: j_decay = " << j_decay << endl ;
      QDP_abort(1) ;
    }

    ReadPath(xml_in, path, Pt_src) ;
    ReadPath(xml_in, path, Sl_src) ;
    ReadPath(xml_in, path, Pt_snk) ;
    ReadPath(xml_in, path, Sl_snk) ;

    ReadPath(xml_in, path, MesonP) ;

    read(xml_in, path + "/num_mom", mom2_max) ;

    // GTF: avg_equiv_mom not part of version 5, true by default
    avg_equiv_mom = true ;

    {
      int input_wvf_kind ;
      read(xml_in, path + "/Wvf_kind", input_wvf_kind) ;
      switch (input_wvf_kind) {
      case 3 :
        Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN ;
        break ;
      default :
        cerr << "Unsupported gauge-invariant Wvf_kind." << endl ;
        cerr << "  Wvf_kind = " << input_wvf_kind << endl ;
        QDP_abort(1) ;
      }
    }

    wvf_param.resize(numKappa) ;
    ReadPath(xml_in, path, wvf_param) ;

    WvfIntPar.resize(numKappa) ;
    ReadPath(xml_in, path, WvfIntPar) ;

    disk_prop.resize(numKappa) ;
    ReadPath(xml_in, path, disk_prop) ;
    for (i=0; i<numKappa; ++i) {
      if (toBool(disk_prop[i] != 1)) {
        cerr << "Only propagators on disk are supported." << endl ;
        cerr << "  disk_prop[" << i << "] = " << disk_prop[i] << endl ;
        QDP_abort(1) ;
      }
    }

    ReadPath(xml_in, path, nrow) ;
    ReadPath(xml_in, path, boundary) ;
    ReadPath(xml_in, path, t_srce) ;

    // Read in the gauge configuration file name
    ReadPath(xml_in, xml_in_root + "/Cfg", cfg_file) ;
    
    break ;

  /**************************************************************************/
  case 6 :
  /**************************************************************************/

    path = xml_in_root + "/param" ; // push into 'param' group

    cerr << "DEBUG" << endl ;

    {
      int input_Nc ;
      read(xml_in, path + "/Nc", input_Nc) ;

      if (input_Nc != Nc) {
        cerr << "Input parameter Nc=" << input_Nc \
          <<  " different from qdp++ value." << endl ;
        QDP_abort(1) ;
      }

      int input_Nd ;
      read(xml_in, path + "/Nd", input_Nd) ;

      if (input_Nd != Nd) {
        cerr << "Input parameter Nd=" << input_Nd \
          << " different from qdp++ value." << endl ;
        QDP_abort(1) ;
      }

      string ferm_type_str ;
      read(xml_in, path + "/FermTypeP", ferm_type_str) ;
      if (ferm_type_str == "WILSON") {
        FermTypeP = FERM_TYPE_WILSON ;
      } else {
        FermTypeP = FERM_TYPE_UNKNOWN ;
      }
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (FermTypeP) {
    case FERM_TYPE_WILSON :

      cout << " SPECTRUM_W: Spectroscopy for Wilson fermions" << endl ;

      ReadPath(xml_in, path, numKappa) ;

      Kappa.resize(numKappa) ;
      ReadPath(xml_in, path, Kappa) ;

      for (i=0; i < numKappa; ++i) {
        if (toBool(Kappa[i] < 0.0)) {
          cerr << "Unreasonable value for Kappa." << endl ;
          cerr << "  Kappa[" << i << "] = " << Kappa[i] << endl ;
          QDP_abort(1) ;
        } else {
          cout << " Spectroscopy Kappa: " << Kappa[i] << endl ;
        }
      }

      break ;

    default :
      cerr << "Fermion type not supported." << endl ;
      if (FermTypeP == FERM_TYPE_UNKNOWN) {
        cerr << "  FermTypeP = UNKNOWN" << endl ;
      }
      QDP_abort(1) ;
    }

    {
      string cfg_type_str ;
      read(xml_in, path + "/cfg_type", cfg_type_str) ;
      if (cfg_type_str == "SZIN") {
        cfg_type = CFG_TYPE_SZIN ;
      } else {
        cfg_type = CFG_TYPE_UNKNOWN ;
      }
    }

    ReadPath(xml_in, path, j_decay) ;
    if (j_decay < 0 || j_decay >= Nd) {
      cerr << "Bad value: j_decay = " << j_decay << endl ;
      QDP_abort(1) ;
    }

    ReadPath(xml_in, path, Pt_src) ;
    ReadPath(xml_in, path, Sl_src) ;
    ReadPath(xml_in, path, Pt_snk) ;
    ReadPath(xml_in, path, Sl_snk) ;

    ReadPath(xml_in, path, MesonP) ;

    ReadPath(xml_in, path, mom2_max) ;
    ReadPath(xml_in, path, avg_equiv_mom) ;

    {
      string wvf_kind_str ;
      read(xml_in, path + "/Wvf_kind", wvf_kind_str) ;
      if (wvf_kind_str == "GAUGE_INV_GAUSSIAN") {
        Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN ;
      } else {
        cerr << "Unsupported gauge-invariant Wvf_kind." << endl ;
        cerr << "  Wvf_kind = " << wvf_kind_str << endl ;
        QDP_abort(1) ;
      }
    }

    wvf_param.resize(numKappa) ;
    ReadPath(xml_in, path, wvf_param) ;

    WvfIntPar.resize(numKappa) ;
    ReadPath(xml_in, path, WvfIntPar) ;

    ReadPath(xml_in, path, nrow) ;
    ReadPath(xml_in, path, boundary) ;
    ReadPath(xml_in, path, t_srce) ;

    // Read in the gauge configuration file name
    ReadPath(xml_in, xml_in_root + "/Cfg", cfg_file) ;

    break ;

  /**************************************************************************/
  default :
  /**************************************************************************/

    cerr << "Input parameter version " << version << " unsupported." << endl ;
    QDP_abort(1) ;
  }

  // GTF: ALL NAMELIST INPUT COMPLETED
  xml_in.close() ;

  // Check for valid parameter values
  for (i=0; i<Nd; ++i) {
    if ((nrow[i] % 2) != 0) {
      cerr << "Lattice shape is invalid; odd linear size not allowed." \
        << endl ;
      cerr << "  nrow[" << i << "] = " << nrow[i] << endl ;
      QDP_abort(1) ;
    }
  }

  // Check that only one type of source smearing is specified
  // Must match how the stored propagator was smeared
  if ((Pt_src == true) && (Sl_src == true)) {
    cerr << "Error if Pt_src and Sl_src are both set." << endl ;
    cerr << "Choose the one which matches the stored propagators." << endl ;
    QDP_abort(1) ;
  }

  // Specify lattice size, shape, etc.
  Layout::setLattSize(nrow) ;
  Layout::create() ;

  // Figure out what to do about boundary conditions
  // GTF HACK: only allow periodic boundary conditions
  for (i=0; i<Nd; ++i) {
    if (boundary[i] != 1) {
      cerr << "Only periodic boundary conditions supported." << endl ;
      cerr << "  boundary[" << i << "] = " << boundary[i] << endl ;
      QDP_abort(1) ;
    }
  }

  for (i=0; i<Nd; ++i) {
    if (t_srce[i] < 0 || t_srce[i] >= nrow[i]) {
      cerr << "Quark propagator source coordinate incorrect." << endl ;
      cerr << "t_srce[" << i << "] = " << t_srce[i] << endl ;
      QDP_abort(1) ;
    }
  }

  cout << endl << "     Gauge group: SU(" << Nc << ")" << endl ;

  multi1d<LatticeColorMatrix> u(Nd) ;

  cout << "     volume: " << nrow[0] ;
  for (i=1; i<Nd; ++i) {
    cout << " x " << nrow[i] ;
  }
  cout << endl ;

  // Read in the configuration along with relevant information.
  XMLReader gauge_xml;

  switch (cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, cfg_file);
    break ;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  // Instantiate namelist writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT") ;

  push(xml_out, "spectrum_w");

  // Write out configuration data to namelist output
  push(xml_out, "IO_version") ;
  Write(xml_out, version) ;
  pop(xml_out) ;

  xml_out.flush();

  push(xml_out, "Output_version") ;
  write(xml_out, "out_version", 5) ;
  pop(xml_out) ;

  push(xml_out, "param") ;

  // was "param1"
  switch (FermTypeP) {
  case FERM_TYPE_WILSON :
    write(xml_out, "FermTypeP", "WILSON") ;
    break ;
  default :
    write(xml_out, "FermTypeP", "UNKNOWN") ;
  }
  Write(xml_out, Nd) ;
  Write(xml_out, Nc) ;
  Write(xml_out, Ns) ;
  Write(xml_out, numKappa) ;
  Write(xml_out, Kappa) ;

  // was "param2"

  // was "param3"
  switch (cfg_type) {
  case CFG_TYPE_SZIN :
    write(xml_out, "cfg_type", "SZIN") ;
    break ;
  default :
    write(xml_out, "cfg_type", "UNKNOWN") ;
  }
  Write(xml_out, j_decay) ;

  // was "param4"

  // was "param5"
  Write(xml_out, Pt_src) ;
  Write(xml_out, Sl_src) ;
  Write(xml_out, Pt_snk) ;
  Write(xml_out, Sl_snk) ;

  // was "param6"
  switch (Wvf_kind) {
  case WVF_KIND_GAUGE_INV_GAUSSIAN :
    write(xml_out, "Wvf_kind", "GAUGE_INV_GAUSSIAN") ;
    break ;
  default :
    write(xml_out, "Wvf_kind", "UNKNOWN") ;
  }
  Write(xml_out, wvf_param) ;
  Write(xml_out, WvfIntPar) ;

  // was "param7"
  // GTF: don't bother since disk_prop is always true
  //Write(xml_out, disk_prop) ;

  // was "param8"

  // was "param9"
  Write(xml_out, MesonP) ;

  // was "param10"

  // before seed, write out mom2_max and avg_equiv_mom
  Write(xml_out, mom2_max) ;
  Write(xml_out, avg_equiv_mom) ;

  pop(xml_out) ;

  push(xml_out, "lattis") ;
  Write(xml_out, nrow) ;
  Write(xml_out, boundary) ;
  Write(xml_out, t_srce) ;
  pop(xml_out) ;

  // Write out the config info
  write(xml_out, "config_info", gauge_xml);

  xml_out.flush();

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
  Double w_plaq, s_plaq, t_plaq, link ;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link) ;

  push(xml_out, "Observables") ;
  Write(xml_out, w_plaq) ;
  Write(xml_out, s_plaq) ;
  Write(xml_out, t_plaq) ;
  Write(xml_out, link) ;
  pop(xml_out) ;

  xml_out.flush();

  // Next check the gauge field configuration by reunitarizing.
  {
    LatticeBoolean lbad ;
    int numbad ;
    multi1d<LatticeColorMatrix> u_tmp(Nd) ;
    u_tmp = u ;
    int mu ;
    for (mu=0; mu < Nd; ++mu) {
      reunit(u_tmp[mu], lbad, numbad, REUNITARIZE_ERROR) ;
    }
  }

  // Initialize the slow Fourier transform phases
  SftMom phases(mom2_max, avg_equiv_mom, j_decay) ;

  // Keep an array of all the xml output buffers
  XMLArrayWriter xml_array(xml_out,numKappa);
  push(xml_array, "Wilson_hadron_measurements");

  // Now loop over the various fermion masses
  for (int loop=0; loop < numKappa; ++loop)
  {
    // Read the quark propagator
    LatticePropagator quark_propagator ;
    {
      XMLReader prop_xml;
      stringstream prop_file ;
      prop_file << "propagator_" << loop ;
      readSzinQprop(prop_xml, quark_propagator, prop_file.str()) ;
    }

    push(xml_array);         // next array element - name auto-written
    Write(xml_array, loop) ;
    write(xml_array, "Kappa_mes", Kappa[loop]) ;
    Write(xml_array, t_srce) ;

    // Do Point sink hadrons and the currents first
    if (MesonP == true) {
      // Construct {Point|Shell}-Point mesons, if desired
      if (Pt_snk == true) {
        if (Pt_src == true) {
          mesons(quark_propagator, quark_propagator, phases, t_srce[j_decay],
                 xml_array, "Point_Point_Wilson_Mesons") ;
        } else if (Sl_src == true) {
          mesons(quark_propagator, quark_propagator, phases, t_srce[j_decay],
                 xml_array, "Shell_Point_Wilson_Mesons") ;
        }
      } // end if (Pt_snk)

      // Convolute the quark propagator with the sink smearing function.
      // Make a copy of the quark propagator and then overwrite it with
      // the convolution. 
      if (Sl_snk == true) {
        LatticePropagator quark_prop_smr ;
        quark_prop_smr = quark_propagator ;
        sink_smear2(u, quark_prop_smr, Wvf_kind, wvf_param[loop],
                    WvfIntPar[loop], j_decay) ;
        if (Pt_src == true) {
          mesons(quark_prop_smr, quark_prop_smr, phases, t_srce[j_decay],
                 xml_array, "Point_Shell_Wilson_Mesons") ;
        } else if (Sl_src == true) {
          mesons(quark_prop_smr, quark_prop_smr, phases, t_srce[j_decay],
                 xml_array, "Shell_Shell_Wilson_Mesons") ;
        }
      } // end if (Sl_snk)

    } // end if (MesonP)

    pop(xml_array) ;  // array element

  } // end for(loop)

  pop(xml_array) ;  // Wilson_spectroscopy
  pop(xml_out) ;  // spectrum_w

  xml_out.close() ;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
