// $Id: spectrum_w.cc,v 1.47 2005-02-28 03:34:47 edwards Exp $
/*! \file
 * \brief Main code for spectrum measurements
 */

#include "chroma.h"

using namespace Chroma;


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

  bool HybMesP;            // Hybrid meson spectroscopy
  int  numb_sm;            // number of smearing levels for E- and B-fields
  Real fact_sm;            // Smearing factor for "smeared" E- and B-fields 

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

  param.HybMesP = false;

  switch (version) 
  {
  case 9:
    param.Wl_snk = false;
    break;

  case 10:
    read(paramtop, "Wl_snk", param.Wl_snk);
    break;

  case 11:
    read(paramtop, "Wl_snk", param.Wl_snk);
    read(paramtop, "HybMesP", param.HybMesP);
    break;

  default:
    QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
    QDP_abort(1);
  }

  if (paramtop.count("HybMesP") != 0)
  {
    // Must read f_mu-nu smearing params
    read(paramtop, "fact_sm", param.fact_sm);
    read(paramtop, "numb_sm", param.numb_sm);
  }
  else
  {
    param.fact_sm = zero;
    param.numb_sm = 0;
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
  ChromaInitialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  Spectrum_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("./DATA");

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
  XMLFileWriter& xml_out = TheXMLOutputWriter::Instance();
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
  MesPlq(xml_out, "Observables", u);
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
    multi1d<int> boundary = getFermActBoundary(prop_header.fermact);
    multi1d<int> t_source = source_header.t_source;

    // Hunt around to find the mass
    // NOTE: this may be problematic in the future if actions are used with no
    // clear def. of a Mass
    std::istringstream  xml_s(prop_header.fermact);
    XMLReader  fermacttop(xml_s);
    const string fermact_path = "/FermionAction";
    string fermact;
    Real Mass;

    QDPIO::cout << "Try action and mass" << endl;
    try
    {
      XMLReader top(fermacttop, fermact_path);

      read(top, "FermAct", fermact);

      // Yuk - need to hop some hoops. This should be isolated.
      if (top.count("Mass") != 0) 
      {
        read(top, "Mass", Mass);
      }
      else if (top.count("Kappa") != 0)
      {
        Real Kappa;
        read(top, "Kappa", Kappa);
        Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
      }
      else if (top.count("m_q") != 0) 
      {
        read(top, "m_q", Mass);
      }
      else
      {
        QDPIO::cerr << "Neither Mass nor Kappa found" << endl;
        throw std::string("Neither Mass nor Kappa found");
      }
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error reading fermact or mass: " << e << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << "FermAct = " << fermact << endl;
    QDPIO::cout << "Mass = " << Mass << endl;

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


    // Next do the hybrid mesons
    if (input.param.HybMesP) 
    {
      /* Smear the gauge fields to construct smeared E- and B-fields */
      int BlkMax = 100;
      Real BlkAccu = fuzz;
      BlkAccu *= 0.01;
      multi1d<LatticeColorMatrix> f;
      multi1d<LatticeColorMatrix> u_smr = u;
      multi1d<LatticeColorMatrix> u_tmp(Nd);

      for(int i=0; i < input.param.numb_sm; ++i)
      {
	for(int mu = 0; mu < Nd; ++mu)
	{
	  // Smear all directions, including time
	  APE_Smear(u_smr, u_tmp[mu], mu, 0, 
		    input.param.fact_sm, BlkAccu, BlkMax, 
		    j_decay);
	}

	u_smr = u_tmp;
      }

      mesField(f, u_smr);  // compute F_munu fields

      // Make traceless (is already anti-hermitian)
      for(int i = 0; i < f.size(); ++i) {
	taproj(f[i]);
      }

      // Construct {Point|Shell}-Point hybrid mesons, if desired
      if (input.param.Pt_snk) 
      {
	if (Pt_src)
	  hybmeson(f, u_smr, quark_propagator, quark_propagator, phases, t_source,
		   xml_array, "Point_Point_Wilson_Hybrid_Mesons");
        
	if (Sl_src)
	  hybmeson(f, u_smr, quark_propagator, quark_propagator, phases, t_source,
		   xml_array, "Shell_Point_Wilson_Hybrid_Mesons");
        
	if (Wl_src)
	  QDP_error_exit("Wall-source hybrid mesons not supported");

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
	  hybmeson(f, u_smr, quark_prop_smr, quark_prop_smr, phases, t_source,
		   xml_array, "Point_Shell_Wilson_Hybrid_Mesons");

	if (Sl_src)
	  hybmeson(f, u_smr, quark_prop_smr, quark_prop_smr, phases, t_source,
		   xml_array, "Shell_Shell_Wilson_Hybrid_Mesons");

	if (Wl_src)
	  QDP_error_exit("Wall-source hybrid mesons not supported");
      } // end if (Sl_snk)

      // Wall sink
      if (input.param.Wl_snk) 
      {
	QDP_error_exit("Wall-sink not supported in hybrid mesons");
      } // end if (Wl_snk)

    } // end if (HybMesP)



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

  END_CODE();

  // Time to bolt
  ChromaFinalize();

  exit(0);
}
