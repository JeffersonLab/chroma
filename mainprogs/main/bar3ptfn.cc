// $Id: bar3ptfn.cc,v 1.9 2003-05-30 02:37:40 flemingg Exp $
//
// $Log: bar3ptfn.cc,v $
// Revision 1.9  2003-05-30 02:37:40  flemingg
// A message printed to stdout was printing the wrong thing cut and pasted
// from spectrum_w.cc.  Replaced with the intended thing from szin bar3ptfn.m
//
// Revision 1.8  2003/05/14 06:07:03  flemingg
// Should be done playing around with the input and output formats
// for bar3ptfn with this version.
//
//

#include "chroma.h"

using namespace QDP ;

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

int
main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv) ;

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

  int t_sink ;

  multi1d<int> sink_mom(Nd-1) ;

  int mom2_max ;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.

  WvfKind Wvf_kind ;        // Wave function kind: gauge invariant
  multi1d<Real> wvf_param ; // Array of width's or other parameters
                            //   for "shell" source/sink wave function
  multi1d<int> WvfIntPar ;  // Array of iter numbers to approx. Gaussian or
                            //   terminate CG inversion for Wuppertal smearing

  int numSeq_src ;
  multi1d<int> Seq_src ;

  multi1d<int> nrow(Nd) ;
  multi1d<int> boundary(Nd) ;
  multi1d<int> t_srce(Nd) ;

  string cfg_file ;

  Seed seed ; // Random number seed (see SETRN for meaning)

  // Instantiate namelist reader for DATA
  NmlReader nml_in("DATA") ;

  // Read in the IO_version
  push(nml_in, "IO_version") ;
  Read(nml_in, version) ;
  pop(nml_in) ;

  switch (version) {

  /**************************************************************************/
  case 3 :
  /**************************************************************************/

    push(nml_in, "param") ; // push into 'param' group

    {
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

      int ferm_type_int ;
      read(nml_in, "FermTypeP", ferm_type_int) ;
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

      cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl ;

      Read(nml_in, numKappa) ;

      Kappa.resize(numKappa) ;
      Read(nml_in, Kappa) ;

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
      read(nml_in, "cfg_type", input_cfg_type) ;
      switch (input_cfg_type) {
      case 1 :
        cfg_type = CFG_TYPE_SZIN ;
        break ;
      default :
        cfg_type = CFG_TYPE_UNKNOWN ;
      }
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

    Read(nml_in, t_sink) ;

    Read(nml_in, sink_mom) ;

    {
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
    }

    wvf_param.resize(numKappa) ;
    Read(nml_in, wvf_param) ;

    WvfIntPar.resize(numKappa) ;
    Read(nml_in, WvfIntPar) ;

    // Now we read in the information associated with the sequential sources
    Read(nml_in, numSeq_src) ;

    // Now read in the particular Sequential Sources we are evaluating
    Seq_src.resize(numSeq_src) ;
    Read(nml_in, Seq_src) ;

    {
      int seq_src_ctr ;
      for (seq_src_ctr=0; seq_src_ctr<numSeq_src; ++seq_src_ctr) {
        cout << "Computing sequential source of type "
          << Seq_src[seq_src_ctr] << endl ;
      }
    }

    Read(nml_in, nrow) ;

    Read(nml_in, boundary) ;

    Read(nml_in, t_srce) ;

    // default value in SZIN
    mom2_max = 7 ;

    pop(nml_in) ;

    // Read in the gauge configuration file name
    push(nml_in, "Cfg") ;
    Read(nml_in, cfg_file) ;
    pop(nml_in) ;

    break ;

  /**************************************************************************/
  case 4 :
  /**************************************************************************/

    push(nml_in, "param") ; // push into 'param' group

    {
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

      string ferm_type_str ;
      read(nml_in, "FermTypeP", ferm_type_str) ;
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

      cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl ;

      Read(nml_in, numKappa) ;

      Kappa.resize(numKappa) ;
      Read(nml_in, Kappa) ;

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
      read(nml_in, "cfg_type", cfg_type_str) ;
      if (cfg_type_str == "SZIN") {
        cfg_type = CFG_TYPE_SZIN ;
      } else {
        cfg_type = CFG_TYPE_UNKNOWN ;
      }
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

    Read(nml_in, t_sink) ;

    Read(nml_in, sink_mom) ;

    {
      string wvf_kind_str ;
      read(nml_in, "Wvf_kind", wvf_kind_str) ;
      if (wvf_kind_str == "GAUGE_INV_GAUSSIAN") {
        Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN ;
      } else {
        cerr << "Unsupported gauge-invariant Wvf_kind." << endl ;
        cerr << "  Wvf_kind = " << wvf_kind_str << endl ;
        QDP_abort(1) ;
      }
    }

    wvf_param.resize(numKappa) ;
    Read(nml_in, wvf_param) ;

    WvfIntPar.resize(numKappa) ;
    Read(nml_in, WvfIntPar) ;

    // Now we read in the information associated with the sequential sources
    Read(nml_in, numSeq_src) ;

    // Now read in the particular Sequential Sources we are evaluating
    Seq_src.resize(numSeq_src) ;
    Read(nml_in, Seq_src) ;

    {
      int seq_src_ctr ;
      for (seq_src_ctr=0; seq_src_ctr<numSeq_src; ++seq_src_ctr) {
        cout << "Computing sequential source of type "
          << Seq_src[seq_src_ctr] << endl ;
      }
    }

    Read(nml_in, nrow) ;

    Read(nml_in, boundary) ;

    Read(nml_in, t_srce) ;

    Read(nml_in, mom2_max) ;

    pop(nml_in) ;

    // Read in the gauge configuration file name
    push(nml_in, "Cfg") ;
    Read(nml_in, cfg_file) ;
    pop(nml_in) ;

  break ;

  /**************************************************************************/
  default :
  /**************************************************************************/

    cerr << "Input parameter version " << version << " unsupported." << endl ;
    QDP_abort(1) ;
  }

  // GTF: ALL NAMELIST INPUT COMPLETED
  nml_in.close() ;

  // Specify lattice size, shape, etc.
  Layout::setLattSize(nrow) ;
  Layout::create() ;

  // Check for valid parameter values
  for (i=0; i<Nd; ++i) {
    if ((nrow[i] % 2) != 0) {
      cerr << "Lattice shape is invalid; odd linear size not allowed." \
        << endl ;
      cerr << "  nrow[" << i << "] = " << nrow[i] << endl ;
      QDP_abort(1) ;
    }
  }

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

  if (t_sink < 0 || t_sink >= nrow[j_decay]) {
    cerr << "Sink time coordinate incorrect." << endl ;
    cerr << "t_sink = " << t_sink << endl ;
    QDP_abort(1) ;
  }

  cout << endl << "     Gauge group: SU(" << Nc << ")" << endl ;

  // Check for unnecessary multiple occurances of kappas and/or wvf_params
  if (numKappa > 1) {
    if (Sl_src == true) {
      for (i=1; i < numKappa; ++i) {
        for (j=0; j<i; ++j) {
          if (toBool(Kappa[j] == Kappa[i])
              && toBool(wvf_param[j] == wvf_param[i])) {
            cerr << "Same kappa and wvf_param:" << endl ;
            cerr << "  Kappa["     << i << "] = " << Kappa[i]     << endl ;
            cerr << "  wvf_param[" << i << "] = " << wvf_param[i] << endl ;
            cerr << "  Kappa["     << j << "] = " << Kappa[j]     << endl ;
            cerr << "  wvf_param[" << j << "] = " << wvf_param[j] << endl ;
            QDP_abort(1) ;
          }
        }
      }
    } else {
      for (i=1; i < numKappa; ++i) {
        for (j=0; j<i; ++j) {
          if (toBool(Kappa[j] == Kappa[i])) {
            cerr  << "Same kappa without shell source or sink:" << endl ;
            cerr << "  Kappa["     << i << "] = " << Kappa[i]     << endl ;
            cerr << "  Kappa["     << j << "] = " << Kappa[j]     << endl ;
            QDP_abort(1) ;
          }
        }
      }
    }
  }

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
  write(nml_out, "out_version", 4) ;
  pop(nml_out) ;

  push(nml_out, "param") ;

  switch (FermTypeP) {
  case FERM_TYPE_WILSON :
    write(nml_out, "FermTypeP", "WILSON") ;
    break ;
  default :
    write(nml_out, "FermTypeP", "UNKNOWN") ;
  }
  Write(nml_out, Nd) ;
  Write(nml_out, Nc) ;
  Write(nml_out, Ns) ;
  Write(nml_out, numKappa) ;
  Write(nml_out, Kappa) ;

  switch (cfg_type) {
  case CFG_TYPE_SZIN :
    write(nml_out, "cfg_type", "SZIN") ;
    break ;
  default :
    write(nml_out, "cfg_type", "UNKNOWN") ;
  }
  Write(nml_out, j_decay) ;

  Write(nml_out, Pt_src) ;
  Write(nml_out, Sl_src) ;
  Write(nml_out, Pt_snk) ;
  Write(nml_out, Sl_snk) ;

  Write(nml_out, t_sink) ;
  Write(nml_out, sink_mom) ;

  switch (Wvf_kind) {
  case WVF_KIND_GAUGE_INV_GAUSSIAN :
    write(nml_out, "Wvf_kind", "GAUGE_INV_GAUSSIAN") ;
    break ;
  default :
    write(nml_out, "Wvf_kind", "UNKNOWN") ;
  }
  Write(nml_out, wvf_param) ;
  Write(nml_out, WvfIntPar) ;

  Write(nml_out, numSeq_src) ;
  Write(nml_out, Seq_src) ;

  Write(nml_out, mom2_max) ;

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
  multi1d<LatticeColorMatrix> u_tmp(Nd) ;
  u_tmp = u ;
  LatticeBoolean lbad = true ;
  int numbad ;
  int mu ;
  for (mu=0; mu < Nd; ++mu) {
    reunit(u_tmp[mu], lbad, numbad, REUNITARIZE_ERROR) ;
  }

  // Now loop over the various kappas
  int loop ;
  for (loop=0; loop < numKappa; ++loop) {

    // Read the quark propagator
    LatticePropagator quark_propagator ;
    {
      stringstream prop_file ;
      prop_file << "propagator_" << loop ;
      readSzinQprop(quark_propagator, prop_file.str()) ;
    }

    int seq_src_ctr ;
    for (seq_src_ctr = 0; seq_src_ctr < numSeq_src; ++seq_src_ctr) {

      LatticePropagator seq_quark_prop ;
      int seq_src_value = Seq_src[seq_src_ctr] ;

      // Read the sequential propagator
      {
        stringstream prop_file ;
        prop_file << "seqprop_" << loop << "_" << seq_src_value ;
        readSzinQprop(seq_quark_prop, prop_file.str()) ;
      }

      if ((0 <= seq_src_value) && (seq_src_value <= 9)) {
        push(nml_out, "Wilson_Baryon_3Pt_fn_measurements") ;
        Write(nml_out, seq_src_value) ;
        Write(nml_out, t_srce) ;
        Write(nml_out, t_sink) ;
        Write(nml_out, sink_mom) ;
        pop(nml_out) ;
      } else if ((10 <= seq_src_value) && (seq_src_value <= 20)) {
        push(nml_out, "Wilson_Meson_3Pt_fn_measurements") ;
        Write(nml_out, seq_src_value) ;
        Write(nml_out, t_srce) ;
        Write(nml_out, t_sink) ;
        Write(nml_out, sink_mom) ;
        pop(nml_out) ;
      } else {
        QDP_error_exit("Unknown sequential source type", seq_src_value) ;
      }

      // Construct the two-pt function from the source point to the sink
      // using only the seq. quark prop.
      // Take hermitian conjugate of the seq. prop, multiply on both sides
      // with gamma_5 = Gamma(G5) and take the trace
      // Use indexing to pull out precisely the source point.
      int G5 = Ns*Ns-1 ;

      // Contract the sequential propagator with itself
      // to form the 2-pt function at the source.
      // Do "source" smearing, if needed
      LatticePropagator seq_quark_prop_tmp = seq_quark_prop ;

      if (Sl_src == true) {
        sink_smear2(u, seq_quark_prop_tmp, Wvf_kind, wvf_param[loop],
                    WvfIntPar[loop], j_decay) ;
      }

//    GTF HACK: Eliminating seq_hadron doesn't work yet, but should work soon.
#if 0
      Complex seq_hadron_0 =
        peekSite(trace(adj(Gamma(G5)*seq_quark_prop_tmp*Gamma(G5))), t_srce) ;
#else
      LatticeComplex seq_hadron = \
        trace(adj(Gamma(G5)*seq_quark_prop_tmp*Gamma(G5))) ;

      Complex seq_hadron_0 = peekSite(seq_hadron, t_srce) ;
#endif

      push(nml_out,"Wilson_hadron_2Pt_fn") ;
      Write(nml_out, t_srce) ;
      Write(nml_out, t_sink) ;
      Write(nml_out, sink_mom) ;
      Write(nml_out, seq_hadron_0) ;
      pop(nml_out) ;

      // Now the 3pt contractions
      SftMom phases(mom2_max, sink_mom, false, j_decay) ;
      FormFac(u, quark_propagator, seq_quark_prop, phases, t_srce[j_decay],
              nml_out) ;

    } // end loop over sequential sources
  } // end loop over the kappa value

  push(nml_out, "End_Wilson_3Pt_fn_measurements") ;
  Write(nml_out, numSeq_src) ;
  pop(nml_out) ;

  // Close the namelist output file NMLDAT
  nml_out.close() ;

  // Time to bolt
  QDP_finalize() ;

  return 0 ;
}
