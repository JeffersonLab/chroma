// $Id: spectrum_w.cc,v 1.38 2004-07-28 03:08:05 edwards Exp $
//
//! \file
//  \brief Main code for propagator generation
//
//  $Log: spectrum_w.cc,v $
//  Revision 1.38  2004-07-28 03:08:05  edwards
//  Added START/END_CODE to all routines. Changed some to not pass an
//  argument.
//
//  Revision 1.37  2004/06/08 14:19:44  edwards
//  Changed out_version to 11.
//
//  Revision 1.36  2004/04/29 00:55:06  edwards
//  Split tests for source type into a switch statement - allows
//  error detection.
//
//  Revision 1.35  2004/04/28 14:34:43  edwards
//  Moved gauge initialization to calling gaugeStartup().
//
//  Revision 1.34  2004/04/16 19:27:39  bjoo
//  Fixed spectrum_w and seqprop
//
//  Revision 1.33  2004/04/06 04:20:33  edwards
//  Added SZINQIO support.
//
//  Revision 1.32  2004/04/05 16:26:27  edwards
//  Moved init of sftmom down after prop is read. Eliminated read
//  of j_decay from input.
//
//  Revision 1.31  2004/04/05 04:19:13  edwards
//  Added initial support for Wall sources/sinks.
//
//  Revision 1.30  2004/02/26 16:38:22  edwards
//  Added test/write-out of forward_prop  prop_corr.
//
//  Revision 1.29  2004/02/23 03:13:58  edwards
//  Major overhaul of input/output model! Now using EXCLUSIVELY
//  SciDAC propagator format for propagators. Now, Param part of input
//  files directly matches source/sink/propagator/seqprop headers
//  of propagators. All ``known'' input of a propagator is derived
//  from its header(s) and used for subsequent calculations.
//
//  Revision 1.28  2004/02/13 15:27:14  sbasak
//  The p-wave and d-wave part has been stripped off the qqq_w file.
//
//  Revision 1.27  2004/02/11 12:51:35  bjoo
//  Stripped out Read() and Write()
//
//  Revision 1.26  2004/02/03 20:05:12  edwards
//  Removed passing j_decay into curcor2
//
//  Revision 1.25  2004/01/31 23:22:01  edwards
//  Added proginfo call.
//
//  Revision 1.24  2004/01/31 22:43:00  edwards
//  Put in tests of array sizes.
//
//  Revision 1.23  2004/01/29 16:44:36  edwards
//  Removed reading of Nd, Nc, and numKappa. Changed from Kappa to Mass!
//  Can also read Kappa for back compatibility.
//
//  Revision 1.22  2004/01/16 15:05:09  kostas
//  fixed the second occurrence of the shell sink bug
//
//  Revision 1.21  2004/01/16 15:01:06  kostas
//  Corrected the shell sink bug
//
//  Revision 1.20  2004/01/14 22:24:38  kostas
//  added reading capability for NERSC confs
//
//  Revision 1.19  2004/01/06 04:59:14  edwards
//  Standardized the IO.
//
//  Revision 1.18  2004/01/05 21:48:36  edwards
//  Changed WVF_KIND to WVF_TYPE.
//
//  Revision 1.17  2003/12/17 17:34:36  edwards
//  Removed tests of boundary.
//
//  Revision 1.16  2003/10/30 02:32:16  edwards
//  Changed output format and number of vector currents measured.
//
//  Revision 1.15  2003/10/10 17:50:32  edwards
//  Removed some extraneous debugging.
//
//  Revision 1.14  2003/10/10 17:48:02  edwards
//  Added missing time_rev for this version of input_io. Other
//  small tweaks.
//
//  Revision 1.13  2003/10/09 20:32:37  edwards
//  Changed all cout/cerr to QDPIO::cout/cerr. Change QDP_info calls
//  to use QDPIO::cout.
//
//  Revision 1.12  2003/10/02 01:21:15  edwards
//  Small tweaks. Changed input group to be program dependent.
//
//  Revision 1.11  2003/10/01 20:23:46  edwards
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


/*
 * Input 
 */
// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  bool Pt_snk;             // point sink
  bool Sl_snk;             // shell sink
  bool Wl_snk;             // wall sink

  bool MesonP;             // Meson spectroscopy
  bool CurrentP;           // Meson currents
  bool BaryonP;            // Baryons spectroscopy

  bool time_rev;           // Use time reversal in baryon spectroscopy


  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  bool avg_equiv_mom;      // average over equivalent momenta
  WvfKind       wvf_kind;  // Wave function kind: gauge invariant
  multi1d<Real> wvf_param; // Array of width's or other parameters
  //   for "shell" source/sink wave function
  multi1d<int> wvfIntPar;  // Array of iter numbers to approx. Gaussian or
  //   terminate CG inversion for Wuppertal smearing

  multi1d<int> nrow;
};


//! Propagators
struct Prop_t
{
  multi1d<string> prop_files;  // The files are expected to be in SciDAC format!
};


//! Mega-structure of parameters
struct Spectrum_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_files", input.prop_files);
}


//! Reader for parameters
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
  case 9:
    /**************************************************************************/
    param.Wl_snk = false;
    break;

  case 10:
    /**************************************************************************/
    read(paramtop, "Wl_snk", param.Wl_snk);
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "Pt_snk", param.Pt_snk);
  read(paramtop, "Sl_snk", param.Sl_snk);

  read(paramtop, "MesonP", param.MesonP);
  read(paramtop, "CurrentP", param.CurrentP);
  read(paramtop, "BaryonP", param.BaryonP);
  read(paramtop, "time_rev", param.time_rev);

  read(paramtop, "mom2_max", param.mom2_max);
  read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);

  read(paramtop, "wvf_kind", param.wvf_kind);
  read(paramtop, "wvf_param", param.wvf_param);
  read(paramtop, "wvfIntPar", param.wvfIntPar);

  if (param.wvf_param.size() != param.wvfIntPar.size())
  {
    QDPIO::cerr << "wvf_param size inconsistent with wvfintpar size" << endl;
    QDP_abort(1);
  }

  read(paramtop, "nrow", param.nrow);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, Spectrum_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the propagator(s) info
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading prop data: " << e << endl;
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

  START_CODE();

  // Input parameter structure
  Spectrum_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/spectrum_w", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " SPECTRUM_W: Spectroscopy for Wilson fermions" << endl;

  /*
   * Sanity checks
   */
  if (input.param.wvf_param.size() != input.prop.prop_files.size())
  {
    QDPIO::cerr << "wvf_param size inconsistent with prop_files size" << endl;
    QDP_abort(1);
  }

  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;

  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "spectrum_w");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 11);
  pop(xml_out);

  xml_out.flush();


  // First calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  multi1d<DComplex> pollp(Nd);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u, pollp[mu], mu);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  write(xml_out, "pollp", pollp);
  pop(xml_out);

  xml_out.flush();

  // Keep an array of all the xml output buffers
  XMLArrayWriter xml_array(xml_out,input.prop.prop_files.size());
  push(xml_array, "Wilson_hadron_measurements");


  // Now loop over the various fermion masses
  for (int loop=0; loop < input.prop.prop_files.size(); ++loop)
  {
    // Read the quark propagator and extract headers
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSource_t source_header;
    {
      XMLReader prop_file_xml, prop_record_xml;
      readQprop(prop_file_xml, 
		prop_record_xml, quark_propagator,
		input.prop.prop_files[loop], QDPIO_SERIAL);

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
    int j_decay = source_header.j_decay;
    Real Mass   = prop_header.FermActHandle->getMass();
    multi1d<int> boundary = prop_header.boundary;
    multi1d<int> t_source = source_header.t_source;

    // Flags
    int t0      = t_source[j_decay];
    int bc_spec = boundary[j_decay];

    // Initialize the slow Fourier transform phases
    SftMom phases(input.param.mom2_max, input.param.avg_equiv_mom, j_decay);

    // Keep a copy of the phases with NO momenta
    SftMom phases_nomom(0, true, j_decay);

    // Next array element - name auto-written
    push(xml_array);
    write(xml_array, "loop", loop);
    write(xml_array, "Mass_mes", Mass);
    write(xml_array, "t_source", t_source);

    // Save prop input
    write(xml_array, "ForwardProp", prop_header);
    write(xml_array, "PropSource", source_header);

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						   phases.getSet());

      push(xml_array, "Forward_prop_correlator");
      write(xml_array, "forward_prop_corr", forward_prop_corr);
      pop(xml_array);
    }

    // Determine what kind of source to use
    bool Pt_src = false;
    bool Sl_src = false;
    bool Wl_src = false;

    switch (source_header.source_type)
    {
    case SRC_TYPE_POINT_SOURCE:
      Pt_src = true;
      break;

    case SRC_TYPE_SHELL_SOURCE:
      Sl_src = true;
      break;

    case SRC_TYPE_WALL_SOURCE:
      Wl_src = true;
      break;

    default:
      QDPIO::cerr << "Unsupported source type" << endl;
      QDP_abort(1);
    }


    // Do the mesons first
    if (input.param.MesonP) 
    {
      // Construct {Point|Shell}-Point mesons, if desired
      if (input.param.Pt_snk) 
      {
	if (Pt_src)
	  mesons(quark_propagator, quark_propagator, phases, t0,
		 xml_array, "Point_Point_Wilson_Mesons");
        
	if (Sl_src)
	  mesons(quark_propagator, quark_propagator, phases, t0,
		 xml_array, "Shell_Point_Wilson_Mesons");
        
	if (Wl_src)
	  mesons(quark_propagator, quark_propagator, phases_nomom, t0,
		 xml_array, "Wall_Point_Wilson_Mesons");
      } // end if (Pt_snk)

      // Convolute the quark propagator with the sink smearing function.
      // Make a copy of the quark propagator and then overwrite it with
      // the convolution. 
      if (input.param.Sl_snk) 
      {
	LatticePropagator quark_prop_smr;
	quark_prop_smr = quark_propagator;
	sink_smear2(u, quark_prop_smr, 
		    input.param.wvf_kind, 
		    input.param.wvf_param[loop],
		    input.param.wvfIntPar[loop], 
		    j_decay);

	if (Pt_src)
	  mesons(quark_prop_smr, quark_prop_smr, phases, t0,
		 xml_array, "Point_Shell_Wilson_Mesons");

	if (Sl_src)
	  mesons(quark_prop_smr, quark_prop_smr, phases, t0,
		 xml_array, "Shell_Shell_Wilson_Mesons");

	if (Wl_src)
	  mesons(quark_prop_smr, quark_prop_smr, phases_nomom, t0,
		 xml_array, "Wall_Shell_Wilson_Mesons");
      } // end if (Sl_snk)

      // Wall sink
      if (input.param.Wl_snk) 
      {
	LatticePropagator wall_quark_prop;
	wall_qprop(wall_quark_prop, quark_propagator, phases_nomom);

	if (Pt_src)
	  mesons(wall_quark_prop, wall_quark_prop, phases_nomom, t0,
		 xml_array, "Point_Wall_Wilson_Mesons");

	if (Sl_src)
	  mesons(wall_quark_prop, wall_quark_prop, phases_nomom, t0,
		 xml_array, "Shell_Wall_Wilson_Mesons");

	if (Wl_src)
	  mesons(wall_quark_prop, wall_quark_prop, phases_nomom, t0,
		 xml_array, "Wall_Wall_Wilson_Mesons");
      } // end if (Wl_snk)

    } // end if (MesonP)


    // Do the currents next
    if (input.param.CurrentP) 
    {
      // Construct the rho vector-current and the pion axial current divergence
      if (Pt_src)
	curcor2(u, quark_propagator, quark_propagator, phases, 
		t0, 3,
		xml_array, "Point_Point_Meson_Currents");
        
      if (Sl_src)
	curcor2(u, quark_propagator, quark_propagator, phases, 
		t0, 3,
		xml_array, "Shell_Point_Meson_Currents");
        
      if (Wl_src)
	curcor2(u, quark_propagator, quark_propagator, phases_nomom, 
		t0, 3,
		xml_array, "Wall_Point_Meson_Currents");
    } // end if (CurrentP)


    // Do the baryons
    if (input.param.BaryonP) 
    {
      // Construct {Point|Shell}-Point mesons, if desired
      if (input.param.Pt_snk) 
      {
	if (Pt_src)
	  baryon(quark_propagator, phases, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Point_Point_Wilson_Baryons");
        
	if (Sl_src)
	  baryon(quark_propagator, phases, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Shell_Point_Wilson_Baryons");
        
	if (Wl_src)
	  baryon(quark_propagator, phases_nomom, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Wall_Point_Wilson_Baryons");
      } // end if (Pt_snk)

      // Convolute the quark propagator with the sink smearing function.
      // Make a copy of the quark propagator and then overwrite it with
      // the convolution. 
      if (input.param.Sl_snk) 
      {
	LatticePropagator quark_prop_smr;
	quark_prop_smr = quark_propagator;
	sink_smear2(u, quark_prop_smr, 
		    input.param.wvf_kind, 
		    input.param.wvf_param[loop],
		    input.param.wvfIntPar[loop], 
		    j_decay);
	if (Pt_src)
	  baryon(quark_prop_smr, phases, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Point_Shell_Wilson_Baryons");
        
	if (Sl_src)
	  baryon(quark_prop_smr, phases, 
		 t0, bc_spec, input.param.time_rev, 	
		 xml_array, "Shell_Shell_Wilson_Baryons");
        
	if (Wl_src)
	  baryon(quark_prop_smr, phases_nomom, 
		 t0, bc_spec, input.param.time_rev, 	
		 xml_array, "Wall_Shell_Wilson_Baryons");
      } // end if (Sl_snk)

      // Wall sink
      if (input.param.Wl_snk) 
      {
	LatticePropagator wall_quark_prop;
	wall_qprop(wall_quark_prop, quark_propagator, phases_nomom);

	if (Pt_src)
	  baryon(wall_quark_prop, phases_nomom, 
		 t0, bc_spec, input.param.time_rev, 
		 xml_array, "Point_Wall_Wilson_Baryons");
        
	if (Sl_src)
	  baryon(wall_quark_prop, phases_nomom, 
		 t0, bc_spec, input.param.time_rev, 	
		 xml_array, "Shell_Wall_Wilson_Baryons");
        
	if (Wl_src)
	  baryon(wall_quark_prop, phases_nomom, 
		 t0, bc_spec, input.param.time_rev, 	
		 xml_array, "Wall_Wall_Wilson_Baryons");
      } // end if (Wl_snk)

    } // end if (BaryonP)

    pop(xml_array);  // array element

  } // end for(loop)

  pop(xml_array);  // Wilson_spectroscopy
  pop(xml_out);  // spectrum_w

  xml_out.close();
  xml_in.close();

  END_CODE();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
