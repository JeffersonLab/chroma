// $Id: bar3ptfn.cc,v 1.2 2003-04-26 06:09:04 flemingg Exp $
//
// $Log: bar3ptfn.cc,v $
// Revision 1.2  2003-04-26 06:09:04  flemingg
// Now uses NmlReader support for bool's and string's.  Main limitations are
// that it still only allows periodic boundary conditions for fermions and
// that mom2_max is hard-coded to be 7, sort of matching szin.  Probably
// should add something to be read from namelist input to allow a different
// value of mom2_max.  Also, it might be nice to decide if we really need
// to read in and write out a whole bunch of namelist key-value pairs
// that are not used for anything else.  Now might be a good time to clean
// up some accidents of history.
//
//

#include "chroma.h"
#include "qdp_util.h" // for readSzin()

using namespace QDP ;

int
main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv) ;

  int i, j ;

  Seed seed ; // Random number seed (see SETRN for meaning)

  // Instantiate namelist reader for DATA
  NmlReader nml_in("DATA") ;

  // Instantiate namelist writer for NMLDAT
  NmlWriter nml_out("NMLDAT") ;

  // Read in the IO_version (checked later)
  int version ; // input parameter version
  push(nml_in, "IO_version") ;
  Read(nml_in, version) ;
  pop(nml_in) ;

  // push into param group
  push(nml_in, "param") ;

  int FermTypeP ;
  Read(nml_in, FermTypeP) ;

  int numKappa ; // Number of Kappa's passed to hadron
  Read(nml_in, numKappa) ;

  int input_Nd ;
  read(nml_in, "Nd", input_Nd) ;

  if (input_Nd != Nd) {
    QDP_error_exit("Input parameter Nd different from qdp++ value.") ;
  }

  int input_Nc ;
  read(nml_in, "Nc", input_Nc) ;

  if (input_Nc != Nc) {
    QDP_error_exit("Input parameter Nc different from qdp++ value.") ;
  }

  multi1d<Real> Kappa(numKappa) ; // Array of Kappa's passed to hadron

  switch (FermTypeP) {
  case 1 :             // GTF HACK, change once WILSON_FERMIONS enum defined
    cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl ;

    Read(nml_in, Kappa) ;

    for (i = 0; i < numKappa; ++i) {
      cout << " Spectroscopy Kappa: " << Kappa[i] << endl ;
    }

    push(nml_out, "IO_version") ;
    Write(nml_out, version) ;
    pop(nml_out) ;
    push(nml_out, "Output_version") ;
    write(nml_out, "out_version", 3) ;
    pop(nml_out) ;

    break ;
  default :
    QDP_error_exit("Fermion type incorrect", FermTypeP) ;
  }

  Real          OverMass ;
  int           RatPolyDeg ;
  Real          PolyArgResc ;
  int           NWilsVec ;

  int           cfg_type ;
  int           j_decay ;

  int           Z3_src ; // Z3 random value source at 2^Z3_src pts/dir

  bool          Pt_src ; // Point source
  bool          Sl_src ; // Shell source

  bool          Pt_snk ; // Point sink
  bool          Sl_snk ; // Shell sink
  int           t_sink ;
  multi1d<int>  sink_mom(Nd-1) ;

  int           InvType ;
  int           FermAct ;
  Real          H_parity ;
  Real          ClovCoeff ;
  Real          u0 ;
  Real          MRover ;

  // GTF ???: Why is MaxCG 'int' and not 'multi1d<int> MaxCG(numKappa)'
  int           MaxCG ;
  multi1d<Real> RsdCG(numKappa) ; // CG accuracy

  int           Wvf_kind ; // Wave function kind: gauge invariant

  // Array of width's or other parameters for "shell" source/sink wave function
  multi1d<Real> wvf_param(numKappa) ;

  // Array of iter numbers to approx. Gaussian or terminate CG inversion
  // for Wuppertal smearing
  multi1d<int>  WvfIntPar(numKappa) ;

  int           numSeq_src ; // The total number of sequential sources

  int           numGamma ;

  // An array containing a list of the sequential source
  multi1d<int>  Seq_src ;

  // A list of the gamma matrices in the usual DeGrand-Rossi encoding
  multi1d<int>  Gamma_list ;


  switch (version) {
  case 3 :

    Read(nml_in, OverMass) ;
    Read(nml_in, RatPolyDeg) ;
    Read(nml_in, PolyArgResc) ;
    Read(nml_in, NWilsVec) ;

    Read(nml_in, cfg_type) ; // Configuration type - szin, Illinois staggered
    Read(nml_in, j_decay) ;  // Direction to measure propagators

    Read(nml_in, Z3_src) ;

    Read(nml_in, Pt_src) ;
    Read(nml_in, Sl_src) ;

    Read(nml_in, Pt_snk) ;
    Read(nml_in, Sl_snk) ;
    // Time coordinate of the sequential quark propagator source (ie. the sink)
    Read(nml_in, t_sink) ;
    Read(nml_in, sink_mom) ; // Sink hadron momentum

    Read(nml_in, InvType) ;
    Read(nml_in, FermAct) ;
    Read(nml_in, H_parity) ;
    Read(nml_in, ClovCoeff) ;
    Read(nml_in, u0) ;
    Read(nml_in, MRover) ;

    Read(nml_in, MaxCG) ;
    Read(nml_in, RsdCG) ;
   
    Read(nml_in, Wvf_kind) ;

    Read(nml_in, wvf_param) ;
   
    Read(nml_in, WvfIntPar) ;

    // Now we read in the information associated with the sequential sources
    Read(nml_in, numSeq_src) ;

    // Now the information associated with the gamma matrices
    Read(nml_in, numGamma) ;

    // Now read in the particular Sequential Sources we are evaluating
    Seq_src.resize(numSeq_src) ;
    Read(nml_in, Seq_src) ;

    int seq_src_ctr ;
    for (seq_src_ctr=0; seq_src_ctr<numSeq_src; ++seq_src_ctr) {
      cout << "Computing sequential source of type "
        << Seq_src[seq_src_ctr] << endl ;
    }

    // Read in the list of gamma matrices
    Gamma_list.resize(numGamma) ;
    Read(nml_in, Gamma_list) ;
   
    break ;
  default :
    QDP_error_exit("Unsupported input parameter version", version) ;
  }

  // Specify lattice size, shape, etc.
  multi1d<int> nrow(Nd) ;
  Read(nml_in, nrow) ;
  Layout::setLattSize(nrow) ;
  Layout::create() ;

  multi1d<int> boundary(Nd) ;
  Read(nml_in, boundary) ;

  // Coordinates of the quark propagator source.
  multi1d<int> t_srce(Nd) ;
  Read(nml_in, t_srce) ;

  // pop out of param group
  pop(nml_in) ;

  // Read in the gauge configuration file name
  push(nml_in, "Cfg") ;
  string cfg_file ;
  Read(nml_in, cfg_file) ;
  pop(nml_in) ;

  // Check for valid parameter values
  for (i=0; i<Nd; ++i) {
    if ((nrow[i] % 2) != 0) {
      QDP_error_exit("lattice shape is invalid odd linear size not allowed",
        i, nrow[i]) ;
    }
  }

  // Figure out what to do about boundary conditions
  // GTF HACK: only allow periodic boundary conditions
  for (i=0; i<Nd; ++i) {
    if (boundary[i] != 1) {
      QDP_error_exit("Only periodic boundary conditions supported") ;
    }
  }

  for (i=0; i<Nd; ++i) {
    if (t_srce[i] < 0 || t_srce[i] >= nrow[i]) {
      QDP_error_exit("quark propagator source coordinate incorrect",
        i, t_srce[i]) ;
    }
  }

  if (t_sink < 0 || t_sink >= nrow[j_decay]) {
    QDP_error_exit("sink time coordinate incorrect", t_sink) ;
  }

  cout << endl << "     Gauge group: SU(" << Nc << ")" << endl ;

  for (i=0; i < numKappa; ++i) {
    if (toBool(Kappa[i] < 0.0)) {
      QDP_error_exit("unreasonable value for Kappa", i, toDouble(Kappa[i])) ;
    }
  }

// GTF ???: What to use for SMALL and HALF
//
//for (i=0; i < numKappa; ++i) {
//  if (MaxCG < 0 || RsdCG[i] < SMALL || RsdCG[i] > HALF) {
//    QDP_error_exit("unreasonable CG parameters", MaxCG, RsdCG[i]) ;
//  }
//}

  // Check for unnecessary multiple occurances of kappas and/or wvf_params
  if (numKappa > 1) {
    if (Sl_src == true) {
      for (i=1; i < numKappa; ++i) {
        for (j=0; j<i; ++j) {
          if (toBool(Kappa[j] == Kappa[i])
              && toBool(wvf_param[j] == wvf_param[i])) {
            QDP_error_exit("Same kappa and wvf_param",
              i, toDouble(Kappa[i]), toDouble(wvf_param[i]),
              j, toDouble(Kappa[j]), toDouble(wvf_param[j])) ;
          }
        }
      }
    } else {
      for (i=1; i < numKappa; ++i) {
        for (j=0; j<i; ++j) {
          if (toBool(Kappa[j] == Kappa[i])) {
            QDP_error_exit("Same kappa without shell source or sink",
              i, toDouble(Kappa[i]), j, toDouble(Kappa[j])) ;
          }
        }
      }
    }
  }

  // Initialize neighbour communication and boundary conditions
  // INITIALIZE_GEOMETRY
  // INITIALIZE_BOUNDARY

  multi1d<LatticeColorMatrix> u(Nd) ;

  cout << "     volume: " << nrow[0] ;
  for (i=1; i<Nd; ++i) {
    cout << " x " << nrow[i] ;
  }
  cout << endl ;

  // Read in the configuration along with relevant information.
  switch (cfg_type) {
  case 1 :
    readSzin2(u, (char *)cfg_file.c_str(), seed) ;
    break ;
  default :
    QDP_error_exit("configuration type is unsupported", cfg_type) ;
  }

  // Finished with READ_NAMELIST input
  nml_in.close() ;

  // The call to setrn MUST go after setbc
  RNG::setrn(seed) ;

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

  push(nml_out, "param1") ;
  Write(nml_out, FermTypeP) ;
  Write(nml_out, Nd) ;
  Write(nml_out, Nc) ;
  Write(nml_out, Ns) ;
  Write(nml_out, numKappa) ;
  Write(nml_out, Kappa) ;
  pop(nml_out) ;

  push(nml_out, "param2") ;
  Write(nml_out, FermAct) ;
  Write(nml_out, OverMass) ;
  Write(nml_out, RatPolyDeg) ;
  Write(nml_out, PolyArgResc) ;
  Write(nml_out, NWilsVec) ;
  Write(nml_out, H_parity) ;
  Write(nml_out, ClovCoeff) ;
  write(nml_out, "ClovCoeffR", 0) ;
  write(nml_out, "ClovCoeffT", 0) ;
  Write(nml_out, u0) ;
  pop(nml_out) ;

  push(nml_out, "param3") ;
  Write(nml_out, cfg_type) ;
  Write(nml_out, j_decay) ;
  pop(nml_out) ;

  push(nml_out, "param4") ;
  Write(nml_out, MaxCG) ;
  Write(nml_out, RsdCG) ;
  Write(nml_out, InvType) ;
  Write(nml_out, MRover) ;
  pop(nml_out) ;

  push(nml_out, "param5") ;
  Write(nml_out, Pt_src) ;
  Write(nml_out, Sl_src) ;
  Write(nml_out, Z3_src) ;
  pop(nml_out) ;

  push(nml_out, "param6") ;
  Write(nml_out, Pt_snk) ;
  Write(nml_out, Sl_snk) ;
  Write(nml_out, t_sink) ;
  Write(nml_out, sink_mom) ;
  pop(nml_out) ;

  push(nml_out, "param7") ;
  Write(nml_out, Wvf_kind) ;
  Write(nml_out, wvf_param) ;
  Write(nml_out, WvfIntPar) ;
  pop(nml_out) ;

  push(nml_out, "param8") ;
  write(nml_out, "AnisoP", false) ;
  write(nml_out, "t_dir", Nd-1) ;
  write(nml_out, "xi_0", 1) ;
  write(nml_out, "xiF_0", 1) ;
  pop(nml_out) ;

  push(nml_out, "param9") ;
  Write(nml_out, numSeq_src) ;
  Write(nml_out, Seq_src) ;
  Write(nml_out, numGamma) ;
  Write(nml_out, Gamma_list) ;
  pop(nml_out) ;

  push(nml_out, "param10") ;
  Write(nml_out, seed) ;
  pop(nml_out) ;

  push(nml_out, "lattis") ;
  Write(nml_out, nrow) ;
  Write(nml_out, boundary) ;
  Write(nml_out, t_srce) ;
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

  // If we require a shell wave function sink type, determine it now:
  int wvf_type = 0 ; // GTF ???: should be enum
  if (Sl_snk == true) {
    switch (Wvf_kind) {
    case 3 :
      wvf_type = OPTION_GAUGE_INV_GAUSSIAN_WVF ;
      break ;
    case 4 :
      wvf_type = OPTION_WUPPERTAL_WVF ;
      break ;
    default :
      QDP_error_exit("Unsupported gauge-invariant Wvf_kind[not 3 or 4]",
        Wvf_kind) ;
    }
  }

  // Allocate the source type
  int src_type ; // GTF ???: should be enum
  if (Pt_src == true) {
    src_type = OPTION_POINT_SOURCE ;
    if (Sl_src == true) {
      cout << " Warning: shell source ignored; do point source" << endl ;
    }
  } else if (Sl_src == true) {
    src_type = OPTION_SHELL_SOURCE ;
  } else {
    QDP_error_exit("Must specify point source or shell source") ;
  }

  // Allocate the sink type
  int snk_type ;
  if (Pt_snk == true) {
    snk_type = OPTION_POINT_SINK ;
    if (Sl_snk == true) {
      cout << " Warning: shell sink ignored; do point sink" << endl ;
    }
  } else if (Sl_snk == true) {
    snk_type = OPTION_SHELL_SINK ;
  } else {
    QDP_error_exit("Must specify point sink or shell sink") ;
  }

  // Now loop over the various kappas
  int loop ;
  for (loop=0; loop < numKappa; ++loop) {

    // Read the quark propagator
    LatticePropagator quark_propagator ;
    {
      stringstream prop_file ;
      prop_file << "propagator_" << loop ;
      // GTF HACK: arg 2 should be just string&
      readSzinQprop2(quark_propagator, (char *)prop_file.str().c_str()) ;
    }

    int seq_src_ctr ;
    for (seq_src_ctr = 0; seq_src_ctr < numSeq_src; ++seq_src_ctr) {

      LatticePropagator seq_quark_prop ;
      int seq_src_value = Seq_src[seq_src_ctr] ;

      // Read the sequential propagator
      {
        stringstream prop_file ;
        prop_file << "seqprop_" << loop << "_" << seq_src_value ;
        // GTF HACK: arg 2 should be just string&
        readSzinQprop2(seq_quark_prop, (char *)prop_file.str().c_str()) ;
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
        sink_smear2(u, seq_quark_prop_tmp, wvf_type, wvf_param[loop],
                    WvfIntPar[loop], j_decay) ;
      }

      // GTF HACK: skipping the "Wilson_hadron_2Pt_fn" part for now.

      // Now the 3pt contractions
      // GTF ???: Note that mom2_max = 7 here to match what was done
      // in szin.  Would like to have more momenta. Add something
      // to namelist input file DATA to accomplish this?
      SftMom phases(7, sink_mom, false, j_decay) ;
      FormFac(u, quark_propagator, seq_quark_prop, phases, t_srce[j_decay],
              nml_out) ;

    } // end loop over sequential sources
  } // end loop over the kappa value

  // Close the namelist output file NMLDAT
  nml_out.close() ;

  // Time to bolt
  QDP_finalize() ;

  return 0 ;
}
