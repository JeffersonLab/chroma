// $Id: spectrum_w.cc,v 1.2 2003-05-08 23:02:12 flemingg Exp $
//
//! \file
//  \brief Main code for propagator generation
//
//  $Log: spectrum_w.cc,v $
//  Revision 1.2  2003-05-08 23:02:12  flemingg
//  Initial version of spectrum_w.  It compiles and reproduces szin meson
//  correlation functions up to a sign.  Namelist input and output starting
//  to evolve away from szin model: a work in progress.  Code runs but
//  more testing needed including resolving the cause of the different signs.
//

#include <iostream>
#include <cstdio>

// #define MAIN

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

  Seed seed ;               // Random number seed (see SETRN for meaning)

  // Instantiate namelist reader for DATA
  NmlReader nml_in("DATA") ;

  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  push(nml_in, "IO_version") ;
  Read(nml_in, version) ;
  pop(nml_in) ;

  switch (version) {
  case 5 :
  case 6 :
    push(nml_in, "param") ; // push into 'param' group

    int input_Nc ;
    read(nml_in, "Nc", input_Nc) ;

    if (input_Nc != Nc) {
      cerr << "Input parameter Nc=" << input_Nc \
        <<  " different from qdp++ value." << endl ;
      QDP_abort(1) ;
    }

    int input_Nd ;
    read(nml_in, "Nd", input_Nd) ;

    if (input_Nd != Nd) {
      cerr << "Input parameter Nd=" << input_Nd \
        << " different from qdp++ value." << endl ;
      QDP_abort(1) ;
    }

    pop(nml_in) ; // pop out of 'param' group

    break ;
  default :
    cerr << "Input parameter version " << version << " unsupported." << endl ;
    QDP_abort(1) ;
  }

  // Now determine 'FermTypeP' from namelist input

  switch (version) {
  case 5 :
    push(nml_in, "param") ; // push into param group

    int ferm_type_int ;
    read(nml_in, "FermTypeP", ferm_type_int) ;
    switch (ferm_type_int) {
    case 1 :
      FermTypeP = FERM_TYPE_WILSON ;
      break ;
    default :
      FermTypeP = FERM_TYPE_UNKNOWN ;
    }

    pop(nml_in) ;

    break ;
  case 6 :
    push(nml_in, "param") ; // push into param group

    {
      string ferm_type_str ;
      read(nml_in, "FermTypeP", ferm_type_str) ;
      if (ferm_type_str == "WILSON_FERMIONS") {
        FermTypeP = FERM_TYPE_WILSON ;
      } else {
        FermTypeP = FERM_TYPE_UNKNOWN ;
      }
    }

    pop(nml_in) ;
    break ;
  default :
    cerr << "How did you manage to get here???" << endl ;
    QDP_abort(1) ;
  }

  // Read the fermion masses from the namelist input

  switch (FermTypeP) {
  case FERM_TYPE_WILSON :

    cout << " SPECTRUM_W: Spectroscopy for Wilson fermions" << endl ;

    push(nml_in, "param") ;

    Read(nml_in, numKappa) ;

    Kappa.resize(numKappa) ;
    Read(nml_in, Kappa) ;

    for (i=0; i < numKappa; ++i) {
      if (toBool(Kappa[i] < 0.0)) {
        cerr << "Unreasonable value for Kappa." << endl ;
        cerr << "  Kappa[" << i << "] = " << Kappa[i] << endl ;
      } else {
        cout << " Spectroscopy Kappa: " << Kappa[i] << endl ;
      }
    }

    pop(nml_in) ;

    break ;

  default :
    cerr << "Unsupported fermion type" << endl ;
    QDP_abort(1) ;
  }

  // Read the rest of namelist input
  switch (version) {

  case 5 :

    push(nml_in, "param") ;

    int input_cfg_type ;
    read(nml_in, "cfg_type", input_cfg_type) ;
    switch (input_cfg_type) {
    case 1 :
      cfg_type = CFG_TYPE_SZIN ;
      break ;
    default :
      cfg_type = CFG_TYPE_UNKNOWN ;
    }

    Read(nml_in, j_decay) ;
    if (j_decay < 0 || j_decay >= Nd) {
      cerr << "Bad value: j_decay = " << j_decay << endl ;
      QDP_abort(1) ;
    }

    Read(nml_in, Pt_src) ;
    Read(nml_in, Sl_src) ;
    Read(nml_in, Pt_snk) ;
    Read(nml_in, Sl_snk) ;

    Read(nml_in, MesonP) ;

    read(nml_in, "num_mom", mom2_max) ;

    // GTF: avg_equiv_mom not part of version 5, true by default
    avg_equiv_mom = true ;

    int input_wvf_kind ;
    read(nml_in, "Wvf_kind", input_wvf_kind) ;
    switch (input_wvf_kind) {
    case 3 :
      Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN ;
      break ;
    default :
      cerr << "Unsupported gauge-invariant Wvf_kind." << endl ;
      cerr << "  Wvf_kind = " << input_wvf_kind << endl ;
      QDP_abort(1) ;
    }

    wvf_param.resize(numKappa) ;
    Read(nml_in, wvf_param) ;

    WvfIntPar.resize(numKappa) ;
    Read(nml_in, WvfIntPar) ;

    disk_prop.resize(numKappa) ;
    Read(nml_in, disk_prop) ;
    for (i=0; i<numKappa; ++i) {
      if (toBool(disk_prop[i] != 1)) {
        cerr << "Only propagators on disk are supported." << endl ;
        cerr << "  disk_prop[" << i << "] = " << disk_prop[i] << endl ;
        QDP_abort(1) ;
      }
    }

    Read(nml_in, nrow) ;

    Read(nml_in, boundary) ;

    Read(nml_in, t_srce) ;

    pop(nml_in) ;

    // Read in the gauge configuration file name
    push(nml_in, "Cfg") ;
    Read(nml_in, cfg_file) ;
    pop(nml_in) ;

    break ;

  default :
    cerr << "What are you doing in here??? Get out!!!" << endl ;
    QDP_abort(1) ;
  }

  // GTF: ALL NAMELIST INPUT COMPLETED
  nml_in.close() ;

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
  switch (cfg_type) {
  case CFG_TYPE_SZIN :
    readSzin(u, cfg_file, seed) ;
    break ;
  default :
    cerr << "Configuration type is unsupported." << endl ;
    QDP_abort(1) ;
  }

  // The call to setrn MUST go after setbc
  RNG::setrn(seed) ;

  // Instantiate namelist writer for NMLDAT
  NmlWriter nml_out("NMLDAT") ;

  // Write out configuration data to namelist output
  push(nml_out, "IO_version") ;
  Write(nml_out, version) ;
  pop(nml_out) ;

  push(nml_out, "Output_version") ;
  write(nml_out, "out_version", 5) ;
  pop(nml_out) ;

  push(nml_out, "param") ;

  // was "param1"
  Write(nml_out, FermTypeP) ;
  Write(nml_out, Nd) ;
  Write(nml_out, Nc) ;
  Write(nml_out, Ns) ;
  Write(nml_out, numKappa) ;
  Write(nml_out, Kappa) ;

  // was "param2"

  // was "param3"
  if (cfg_type == CFG_TYPE_SZIN) {
    write(nml_out, "cfg_type", "SZIN") ;
  } else {
    write(nml_out, "cfg_type", "UNKNOWN") ;
  }
  Write(nml_out, j_decay) ;

  // was "param4"

  // was "param5"
  Write(nml_out, Pt_src) ;
  Write(nml_out, Sl_src) ;
  Write(nml_out, Pt_snk) ;
  Write(nml_out, Sl_snk) ;

  // was "param6"
  if (Wvf_kind == WVF_KIND_GAUGE_INV_GAUSSIAN) {
    write(nml_out, "Wvf_kind", "GAUGE_INV_GAUSSIAN") ;
  } else {
    write(nml_out, "Wvf_kind", "UNKNOWN") ;
  }
  Write(nml_out, wvf_param) ;
  Write(nml_out, WvfIntPar) ;

  // was "param7"
  // GTF: don't bother since disk_prop is always true
  //Write(nml_out, disk_prop) ;

  // was "param8"

  // was "param9"
  Write(nml_out, MesonP) ;

  // was "param10"

  // was "param11"
  Write(nml_out, seed) ;

  pop(nml_out) ;

  push(nml_out, "lattis") ;
  Write(nml_out, nrow) ;
  Write(nml_out, boundary) ;
  Write(nml_out, t_srce) ;
  pop(nml_out) ;

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
  Double w_plaq, s_plaq, t_plaq, link ;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link) ;

  push(nml_out, "Observables") ;
  Write(nml_out, w_plaq) ;
  Write(nml_out, s_plaq) ;
  Write(nml_out, t_plaq) ;
  Write(nml_out, link) ;
  pop(nml_out) ;

  // Next check the gauge field configuration by reunitarizing.
  LatticeBoolean lbad = true ;
  int numbad ;
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd) ;
    u_tmp = u ;
    int mu ;
    for (mu=0; mu < Nd; ++mu) {
      reunit(u_tmp[mu], lbad, numbad, REUNITARIZE_ERROR) ;
    }
  }

  // Initialize the slow Fourier transform phases
  SftMom phases(mom2_max, avg_equiv_mom, j_decay) ;

  // Now loop over the various fermion masses
  int loop ;
  for (loop=0; loop < numKappa; ++loop) {

    // Read the quark propagator
    LatticePropagator quark_propagator ;
    {
      stringstream prop_file ;
      prop_file << "propagator_" << loop ;
      readSzinQprop(quark_propagator, prop_file.str()) ;
    }

    push(nml_out, "Wilson_hadron_measurements") ;
    Write(nml_out, loop) ;
    write(nml_out, "Kappa_mes", Kappa[loop]) ;
    Write(nml_out, t_srce) ;
    pop(nml_out) ;

    // Do Point sink hadrons and the currents first
    if (MesonP == true) {
      // Construct {Point|Shell}-Point mesons, if desired
      if (Pt_snk == true) {
        if (Pt_src == true) {
          mesons(quark_propagator, quark_propagator, phases, t_srce[j_decay],
                 nml_out, "Point_Point_Wilson_Mesons") ;
        } else if (Sl_src == true) {
          mesons(quark_propagator, quark_propagator, phases, t_srce[j_decay],
                 nml_out, "Shell_Point_Wilson_Mesons") ;
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
                 nml_out, "Point_Shell_Wilson_Mesons") ;
        } else if (Sl_src == true) {
          mesons(quark_prop_smr, quark_prop_smr, phases, t_srce[j_decay],
                 nml_out, "Shell_Shell_Wilson_Mesons") ;
        }
      } // end if (Sl_snk)

    } // end if (MesonP)

  } // end for(loop)

  nml_out.close() ;

  // Time to bolt
  QDP_finalize();

  return 0;
}
