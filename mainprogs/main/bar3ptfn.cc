// $Id: bar3ptfn.cc,v 1.17 2003-09-11 15:25:32 edwards Exp $
//
// $Log: bar3ptfn.cc,v $
// Revision 1.17  2003-09-11 15:25:32  edwards
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

int
main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  int i, j;

  // Parameters which must be determined from the namelist input
  // and written to the namelist output

  int version;              // input parameter version
  int out_version = 5;       // output version

  FermType FermTypeP;

  // GTF GRIPE: I would prefer masses rather than kappa values here,
  //   but I'll save that for another day.
  int numKappa;            // number of Wilson masses
  multi1d<Real> Kappa;     // array of Wilson mass values

  CfgType cfg_type;        // storage order for stored gauge configuration
  int j_decay;             // direction to measure propagation

  bool Pt_src;             // point source
  bool Sl_src;             // shell source
  bool Pt_snk;             // point sink
  bool Sl_snk;             // shell sink

  int t_sink;

  multi1d<int> sink_mom(Nd-1);

  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.

  WvfKind Wvf_kind;        // Wave function kind: gauge invariant
  multi1d<Real> wvf_param; // Array of width's or other parameters
  //   for "shell" source/sink wave function
  multi1d<int> WvfIntPar;  // Array of iter numbers to approx. Gaussian or
  //   terminate CG inversion for Wuppertal smearing

  int numSeq_src;
  multi1d<int> Seq_src;

  multi1d<int> nrow(Nd);
  multi1d<int> boundary(Nd);
  multi1d<int> t_srce(Nd);

  string cfg_file;

  // Instantiate namelist reader for DATA
  XMLReader xml_in("DATA");
  string xml_in_root = "/bar3ptfn";

  // Read in the IO_version
  read(xml_in, xml_in_root + "/IO_version/version", version);

  switch (version) {

    /**************************************************************************/
  case 3 :
    /**************************************************************************/
    try
    {
      XMLReader xml(xml_in, xml_in_root + "/param"); // push into 'param' group

      {
	int input_Nc;
	read(xml, "Nc", input_Nc);

	if (input_Nc != Nc) {
	  cerr << "Input parameter Nc=" << input_Nc \
	       <<  " different from qdp++ value." << endl;
	  QDP_abort(1);
	}

	int input_Nd;
	read(xml, "Nd", input_Nd);

	if (input_Nd != Nd) {
	  cerr << "Input parameter Nd=" << input_Nd \
	       << " different from qdp++ value." << endl;
	  QDP_abort(1);
	}

	int ferm_type_int;
	read(xml, "FermTypeP", ferm_type_int);
	switch (ferm_type_int) {
	case 1 :
	  FermTypeP = FERM_TYPE_WILSON;
	  break;
	default :
	  FermTypeP = FERM_TYPE_UNKNOWN;
	}
      }

      // GTF NOTE: I'm going to switch on FermTypeP here because I want
      // to leave open the option of treating masses differently.
      switch (FermTypeP) {
      case FERM_TYPE_WILSON :

	cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl;

	Read(xml, numKappa);

	Kappa.resize(numKappa);
	Read(xml, Kappa);

	for (i=0; i < numKappa; ++i) {
	  if (toBool(Kappa[i] < 0.0)) {
	    cerr << "Unreasonable value for Kappa." << endl;
	    cerr << "  Kappa[" << i << "] = " << Kappa[i] << endl;
	    QDP_abort(1);
	  } else {
	    cout << " Spectroscopy Kappa: " << Kappa[i] << endl;
	  }
	}

	break;

      default :
	cerr << "Fermion type not supported." << endl;
	if (FermTypeP == FERM_TYPE_UNKNOWN) {
	  cerr << "  FermTypeP = UNKNOWN" << endl;
	}
	QDP_abort(1);
      }

      {
	int input_cfg_type;
	read(xml, "cfg_type", input_cfg_type);
	switch (input_cfg_type) {
	case 1 :
	  cfg_type = CFG_TYPE_SZIN;
	  break;
	default :
	  cfg_type = CFG_TYPE_UNKNOWN;
	}
      }

      Read(xml, j_decay);
      if (j_decay < 0 || j_decay >= Nd) {
	cerr << "Bad value: j_decay = " << j_decay << endl;
	QDP_abort(1);
      }

      Read(xml, Pt_src);
      Read(xml, Sl_src);
      Read(xml, Pt_snk);
      Read(xml, Sl_snk);

      Read(xml, t_sink);
      Read(xml, sink_mom);

      {
	int input_wvf_kind;
	read(xml, "Wvf_kind", input_wvf_kind);
	switch (input_wvf_kind) {
	case 3 :
	  Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN;
	  break;
	default :
	  cerr << "Unsupported gauge-invariant Wvf_kind." << endl;
	  cerr << "  Wvf_kind = " << input_wvf_kind << endl;
	  QDP_abort(1);
	}
      }

      wvf_param.resize(numKappa);
      Read(xml, wvf_param);

      WvfIntPar.resize(numKappa);
      Read(xml, WvfIntPar);

      // Now we read in the information associated with the sequential sources
      Read(xml, numSeq_src);

      // Now read in the particular Sequential Sources we are evaluating
      Seq_src.resize(numSeq_src);
      Read(xml, Seq_src);

      for (int seq_src_ctr=0; seq_src_ctr<numSeq_src; ++seq_src_ctr) {
        cout << "Computing sequential source of type "
	     << Seq_src[seq_src_ctr] << endl;
      }

      Read(xml, nrow);
      Read(xml, boundary);
      Read(xml, t_srce);

      // default value in SZIN
      mom2_max = 7;

      // Read in the gauge configuration file name
      read(xml_in, xml_in_root + "/Cfg/cfg_file", cfg_file);
    }
    catch (const string& e) 
    {
      cerr << "Error reading DATA: " << e << endl;
      throw;
    }
    break;

    /**************************************************************************/
  case 4 :
    /**************************************************************************/
    try
    {
      XMLReader xml(xml_in, xml_in_root + "/param"); // push into 'param' group

      {
	int input_Nc;
	read(xml, "Nc", input_Nc);

	if (input_Nc != Nc) {
	  cerr << "Input parameter Nc=" << input_Nc \
	       <<  " different from qdp++ value." << endl;
	  QDP_abort(1);
	}

	int input_Nd;
	read(xml, "Nd", input_Nd);

	if (input_Nd != Nd) {
	  cerr << "Input parameter Nd=" << input_Nd \
	       << " different from qdp++ value." << endl;
	  QDP_abort(1);
	}

	string ferm_type_str;
	read(xml, "FermTypeP", ferm_type_str);
	if (ferm_type_str == "WILSON") {
	  FermTypeP = FERM_TYPE_WILSON;
	} else {
	  FermTypeP = FERM_TYPE_UNKNOWN;
	}
      }

      // GTF NOTE: I'm going to switch on FermTypeP here because I want
      // to leave open the option of treating masses differently.
      switch (FermTypeP) {
      case FERM_TYPE_WILSON :

	cout << " FORMFAC: Baryon form factors for Wilson fermions" << endl;

	Read(xml, numKappa);

	Kappa.resize(numKappa);
	Read(xml, Kappa);

	for (i=0; i < numKappa; ++i) {
	  if (toBool(Kappa[i] < 0.0)) {
	    cerr << "Unreasonable value for Kappa." << endl;
	    cerr << "  Kappa[" << i << "] = " << Kappa[i] << endl;
	    QDP_abort(1);
	  } else {
	    cout << " Spectroscopy Kappa: " << Kappa[i] << endl;
	  }
	}

	break;

      default :
	cerr << "Fermion type not supported." << endl;
	if (FermTypeP == FERM_TYPE_UNKNOWN) {
	  cerr << "  FermTypeP = UNKNOWN" << endl;
	}
	QDP_abort(1);
      }

      {
	string cfg_type_str;
	read(xml, "cfg_type", cfg_type_str);
	if (cfg_type_str == "SZIN") {
	  cfg_type = CFG_TYPE_SZIN;
	} else {
	  cfg_type = CFG_TYPE_UNKNOWN;
	}
      }

      Read(xml, j_decay);
      if (j_decay < 0 || j_decay >= Nd) {
	cerr << "Bad value: j_decay = " << j_decay << endl;
	QDP_abort(1);
      }

      Read(xml, Pt_src);
      Read(xml, Sl_src);
      Read(xml, Pt_snk);
      Read(xml, Sl_snk);

      Read(xml, t_sink);
      Read(xml, sink_mom);

      {
	string wvf_kind_str;
	read(xml, "Wvf_kind", wvf_kind_str);
	if (wvf_kind_str == "GAUGE_INV_GAUSSIAN") {
	  Wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN;
	} else {
	  cerr << "Unsupported gauge-invariant Wvf_kind." << endl;
	  cerr << "  Wvf_kind = " << wvf_kind_str << endl;
	  QDP_abort(1);
	}
      }

      wvf_param.resize(numKappa);
      Read(xml, wvf_param);

      WvfIntPar.resize(numKappa);
      Read(xml, WvfIntPar);

      // Now we read in the information associated with the sequential sources
      Read(xml, numSeq_src);

      // Now read in the particular Sequential Sources we are evaluating
      Seq_src.resize(numSeq_src);
      Read(xml, Seq_src);

      for (int seq_src_ctr=0; seq_src_ctr<numSeq_src; ++seq_src_ctr) 
      {
	cout << "Computing sequential source of type "
	     << Seq_src[seq_src_ctr] << endl;
      }

      Read(xml, nrow);
      Read(xml, boundary);
      Read(xml, t_srce);
      Read(xml, mom2_max);

      // Read in the gauge configuration file name
      read(xml_in, xml_in_root + "/Cfg/cfg_file", cfg_file);
    }
    catch (const string& e) 
    {
      cerr << "Error reading DATA: " << e << endl;
      throw;
    }
    break;

    /**************************************************************************/
  default :
    /**************************************************************************/

    cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  xml_in.close();

  // Specify lattice size, shape, etc.
  Layout::setLattSize(nrow);
  Layout::create();

  // Figure out what to do about boundary conditions
  // GTF HACK: only allow periodic boundary conditions
  for (i=0; i<Nd; ++i) {
    if (boundary[i] != 1) {
      cerr << "Only periodic boundary conditions supported." << endl;
      cerr << "  boundary[" << i << "] = " << boundary[i] << endl;
      QDP_abort(1);
    }
  }

  for (i=0; i<Nd; ++i) {
    if (t_srce[i] < 0 || t_srce[i] >= nrow[i]) {
      cerr << "Quark propagator source coordinate incorrect." << endl;
      cerr << "t_srce[" << i << "] = " << t_srce[i] << endl;
      QDP_abort(1);
    }
  }

  if (t_sink < 0 || t_sink >= nrow[j_decay]) {
    cerr << "Sink time coordinate incorrect." << endl;
    cerr << "t_sink = " << t_sink << endl;
    QDP_abort(1);
  }

  cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  // Check for unnecessary multiple occurances of kappas and/or wvf_params
  if (numKappa > 1) {
    if (Sl_src == true) {
      for (i=1; i < numKappa; ++i) {
        for (j=0; j<i; ++j) {
          if (toBool(Kappa[j] == Kappa[i])
              && toBool(wvf_param[j] == wvf_param[i])) {
            cerr << "Same kappa and wvf_param:" << endl;
            cerr << "  Kappa["     << i << "] = " << Kappa[i]     << endl;
            cerr << "  wvf_param[" << i << "] = " << wvf_param[i] << endl;
            cerr << "  Kappa["     << j << "] = " << Kappa[j]     << endl;
            cerr << "  wvf_param[" << j << "] = " << wvf_param[j] << endl;
            QDP_abort(1);
          }
        }
      }
    } else {
      for (i=1; i < numKappa; ++i) {
        for (j=0; j<i; ++j) {
          if (toBool(Kappa[j] == Kappa[i])) {
            cerr  << "Same kappa without shell source or sink:" << endl;
            cerr << "  Kappa["     << i << "] = " << Kappa[i]     << endl;
            cerr << "  Kappa["     << j << "] = " << Kappa[j]     << endl;
            QDP_abort(1);
          }
        }
      }
    }
  }

  multi1d<LatticeColorMatrix> u(Nd);

  cout << "     volume: " << nrow[0];
  for (i=1; i<Nd; ++i) {
    cout << " x " << nrow[i];
  }
  cout << endl;

  // Read in the configuration along with relevant information.
  QDP_info("Attempt to initialize the gauge field");

  XMLReader gauge_xml;
  switch (cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, cfg_file);
    break;
  default :
    cerr << "Configuration type is unsupported." << endl;
    QDP_abort(1);
  }

  QDP_info("Gauge field successfully initialized");

  // Instantiate namelist writer for NMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "bar3ptfn");

  // Write out configuration data to namelist output
  push(xml_out, "IO_version");
  Write(xml_out, version);
  pop(xml_out);

  push(xml_out, "Output_version");
  Write(xml_out, out_version);
  pop(xml_out);

  push(xml_out, "param");

  switch (FermTypeP) {
  case FERM_TYPE_WILSON :
    write(xml_out, "FermTypeP", "WILSON");
    break;
  default :
    write(xml_out, "FermTypeP", "UNKNOWN");
  }
  Write(xml_out, Nd);
  Write(xml_out, Nc);
  Write(xml_out, Ns);
  Write(xml_out, numKappa);
  Write(xml_out, Kappa);

  switch (cfg_type) {
  case CFG_TYPE_SZIN :
    write(xml_out, "cfg_type", "SZIN");
    break;
  default :
    write(xml_out, "cfg_type", "UNKNOWN");
  }
  Write(xml_out, j_decay);

  Write(xml_out, Pt_src);
  Write(xml_out, Sl_src);
  Write(xml_out, Pt_snk);
  Write(xml_out, Sl_snk);

  Write(xml_out, t_sink);
  Write(xml_out, sink_mom);

  switch (Wvf_kind) {
  case WVF_KIND_GAUGE_INV_GAUSSIAN :
    write(xml_out, "Wvf_kind", "GAUGE_INV_GAUSSIAN");
    break;
  default :
    write(xml_out, "Wvf_kind", "UNKNOWN");
  }
  Write(xml_out, wvf_param);
  Write(xml_out, WvfIntPar);

  Write(xml_out, numSeq_src);
  Write(xml_out, Seq_src);

  Write(xml_out, mom2_max);

  pop(xml_out);

  push(xml_out, "lattis");
  Write(xml_out, nrow);
  Write(xml_out, boundary);
  Write(xml_out, t_srce);
  pop(xml_out);

//  xml_out.flush();

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

//  xml_out.flush();

  // Next check the gauge field configuration by reunitarizing.
  multi1d<LatticeColorMatrix> u_tmp(Nd);
  u_tmp = u;
  {
    LatticeBoolean lbad = true;
    int numbad;
    for (int mu=0; mu < Nd; ++mu) 
    {
      reunit(u_tmp[mu], lbad, numbad, REUNITARIZE_ERROR);
    }
  }

  XMLArrayWriter xml_array(xml_out, numKappa);
  push(xml_array, "Wilson_3Pt_fn_measurements");

  // Now loop over the various kappas
  for (int loop=0; loop < numKappa; ++loop) 
  {
    QDP_info("Mass loop = %d",loop);
  
    push(xml_array);
    Write(xml_array, loop);

    // Read the quark propagator
    QDP_info("Attempt to read forward propagator");
  
    LatticePropagator quark_propagator;
    {
      XMLReader prop_xml;
      stringstream prop_file;
      prop_file << "propagator_" << loop;
      readSzinQprop(prop_xml, quark_propagator, prop_file.str());

      write(xml_array, "Forward_prop_info", prop_xml);
    }

    QDP_info("Forward propagator successfully read");
   

    XMLArrayWriter  xml_seq_src(xml_array, numSeq_src);
    push(xml_seq_src, "Sequential_source");

    for (int seq_src_ctr = 0; seq_src_ctr < numSeq_src; ++seq_src_ctr) 
    {
      push(xml_seq_src);
      Write(xml_seq_src, seq_src_ctr);

      LatticePropagator seq_quark_prop;
      int seq_src_value = Seq_src[seq_src_ctr];

      // Read the sequential propagator
      QDP_info("Attempt to read backward propagator");
   
      {
	XMLReader seqprop_xml;
        stringstream prop_file;
        prop_file << "seqprop_" << loop << "_" << seq_src_value;
        readSzinQprop(seqprop_xml, seq_quark_prop, prop_file.str());

	write(xml_array, "Backward_prop_info", seqprop_xml);
      }

      QDP_info("Backward propagator successfully read");
   
      if ((0 <= seq_src_value) && (seq_src_value <= 9)) {
        write(xml_seq_src, "hadron_type", "BARYON");
      } else if ((10 <= seq_src_value) && (seq_src_value <= 20)) {
        write(xml_seq_src, "hadron_type", "MESON");
      } else if ((21 <= seq_src_value) && (seq_src_value <= 30)) {
        write(xml_seq_src, "hadron_type", "BARYON");
      } else {
        QDP_error_exit("Unknown sequential source type", seq_src_value);
      }

      Write(xml_seq_src, seq_src_value);
      Write(xml_seq_src, t_srce);
      Write(xml_seq_src, t_sink);
      Write(xml_seq_src, sink_mom);
	
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

      if (Sl_src == true) {
        sink_smear2(u, seq_quark_prop_tmp, Wvf_kind, wvf_param[loop],
                    WvfIntPar[loop], j_decay);
      }

//    GTF HACK: Eliminating seq_hadron doesn't work yet, but should work soon.
#if 0
      Complex seq_hadron_0 =
        peekSite(trace(adj(Gamma(G5)*seq_quark_prop_tmp*Gamma(G5))), t_srce);
#else
      LatticeComplex seq_hadron = \
        trace(adj(Gamma(G5)*seq_quark_prop_tmp*Gamma(G5)));

      Complex seq_hadron_0 = peekSite(seq_hadron, t_srce);
#endif

      Write(xml_seq_src, seq_hadron_0);

//      xml_seq_src.flush();

      // Now the 3pt contractions
      SftMom phases(mom2_max, sink_mom, false, j_decay);
      FormFac(u, quark_propagator, seq_quark_prop, phases, t_srce[j_decay],
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
  xml_out.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
