// $Id: make_source.cc,v 1.21 2004-02-04 17:41:56 sbasak Exp $
/*! \file
 *  \brief Main code for source generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;


//! Propagator generation
/*! \defgroup propagator Propagator generation
 *  \ingroup main
 *
 * Main program for propagator generation. Here we need some
 * profound and deep discussion of input parameters.
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  multi1d<int> nrow;

  // Read in params
  XMLReader xml_in("DATA");

  int j_decay;
  int version; 		// The input-parameter version

  Real kappa_fake = 3.14159265359;			// Kappa value
  
  WaveStateType wave_state;       // S-wave(0), P-wave(1), D-wave(2)
  int source_direction; // direction

  SourceType source_type;       // POINT_SOURCE, SHELL_SOURCE, WALL_SOURCE
  WvfKind wvf_kind = WVF_KIND_GAUGE_INV_GAUSSIAN;
  			    // shell source kind always GAUGE_INV_GAUSSIAN
  Real wvf_param = 0;		// smearing width
  int WvfIntPar = 0;            // number of iteration for smearing
  int LaplacePower = 0;         // power of laplacian operator: 0, 1, 2, etc.

  Real sm_fact;		// smearing factor
  int sm_numb;		// number of smearing hits
  int BlkMax;		// max iterations in max-ing trace
  Real BlkAccu;		// accuracy of max-ing

  int disp_length;	// displacement length
  int disp_dir;		// displacement direction: x(0),y(1),z(2)

  multi1d<int> t_source(Nd);

  CfgType cfg_type;
  PropType prop_type;

  string source_filename = "";
  string cfg_file;

  string xml_in_root = "/make_source";

  {
    XMLReader inputtop(xml_in, xml_in_root);

    try
    {
      read(inputtop, "IO_version/version", version);
      
      switch(version) 	// The parameters we read in IO version
      {
      case 4:			
      {
	XMLReader paramtop(inputtop, "param");

	read(paramtop, "j_decay", j_decay);
    	read(paramtop, "wave_state", wave_state);       // S-wave, P-wave etc
	read(paramtop, "source_direction", source_direction);

	read(paramtop, "source_type", source_type);     // Point, Gaussian etc

	if (source_type == SRC_TYPE_SHELL_SOURCE)
	{
	  read(paramtop, "wvf_param", wvf_param);
	  read(paramtop, "WvfIntPar", WvfIntPar);
	  read(paramtop, "LaplacePower", LaplacePower);
	  read(paramtop, "sm_fact", sm_fact);
	  read(paramtop, "sm_numb", sm_numb);
	  read(paramtop, "BlkMax", BlkMax);
	  read(paramtop, "BlkAccu", BlkAccu);
	  read(paramtop, "disp_length", disp_length);
	  read(paramtop, "disp_dir", disp_dir);
	}		
	read(paramtop, "t_srce", t_source);

	read(paramtop, "cfg_type", cfg_type);
	read(paramtop, "prop_type", prop_type);

	// Now get the lattice sizes etc
	read(paramtop, "nrow", nrow);

	if (paramtop.count("source_file") != 0)
	  read(paramtop, "source_file", source_filename);
      }
      break;

      default:
	QDP_error_exit("Unknown io version", version);

      }

      // Read in the gauge configuration file name
      read(inputtop, "Cfg/cfg_file", cfg_file);
    }
    catch(const string& e)
    {
      QDP_error_exit("Error reading in make_source: %s", e.c_str());
    }
  }

  Layout::setLattSize(nrow);
  Layout::create();

  switch(source_type){
  case SRC_TYPE_POINT_SOURCE:
    QDPIO::cout << "Point source" << endl;
    break;
  case SRC_TYPE_WALL_SOURCE:
    QDPIO::cout << "Wall source" << endl;
    break;
  case SRC_TYPE_SHELL_SOURCE:
    QDPIO::cout << "Smeared source wvf_param= " << wvf_param <<": WvfIntPar= "
		<< WvfIntPar << endl
		<< "Power of Laplacian operator= " << LaplacePower << endl
		<< "Displacement length= " << disp_length
		<<": Displacement direction= " << disp_dir << endl;
    break;
  default:
    QDPIO::cout << "Unknown source_type" << endl;
    QDP_abort(1);
  }

  // Read a gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;
  Seed seed_old;

  switch (cfg_type) 
  {
  case CFG_TYPE_SZIN:
    readSzin(gauge_xml, u, cfg_file);
    read(gauge_xml, "/szin/seed", seed_old);
    break;

  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, cfg_file);
    seed_old = 11;
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  multi1d<LatticeColorMatrix> u_smr(Nd);
  u_smr = u;
  for(int i=0; i < sm_numb; ++i)
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd);

    for(int mu = 0; mu < Nd; ++mu)
      if ( mu != j_decay )
        APE_Smear(u_smr, u_tmp[mu], mu, 0, sm_fact, BlkAccu, BlkMax, j_decay);
      else
        u_tmp[mu] = u_smr[mu];

    u_smr = u_tmp;
  }



  // Set the rng seed
  QDP::RNG::setrn (seed_old);


  // Useful info

  string xml_filename;

  if (source_filename == "")
  {
    // Set source filename and xml_filename
    switch(wave_state)
    {
  case WAVE_TYPE_S_WAVE:
    xml_filename = "source.xml";
    source_filename = "source_0";
    break;
  case WAVE_TYPE_P_WAVE:
    xml_filename = "p_source.xml";
    source_filename = "p_source_0";
    break;
  case WAVE_TYPE_D_WAVE:    /* added */
    if (source_direction == 12) {
      xml_filename = "dydz_propagator.xml";
      source_filename = "dydz_source_0";
    }
    if (source_direction == 22) {
      xml_filename = "dzdz_propagator.xml";
      source_filename = "dzdz_source_0";
    }
    break;
  default: 
    cerr<<"invaid wave_state\n";
    break;
   } 
  }
  else
  {
    xml_filename = "make_source.xml";
  }


  XMLFileWriter xml_out(xml_filename);
  push(xml_out,"make_source");

  xml_out << xml_in;  // save a copy of the input

  // Write out the config info
  try {
    write(xml_out, "config_info", gauge_xml);
  } 
  catch(const string& e)
  {
    QDP_error_exit("Error writing gauge_xml: %s",e.c_str());
  }

  xml_out.flush();


  //
  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //
  // For this calculation, a smeared source is used. A point
  // source is first constructed and then smeared. If a user
  // only wanted a point source, then remove the smearing stuff
  //
  LatticePropagator quark_source;

  PropHead header;		// Header information
  header.kappa = kappa_fake;
  header.source_smearingparam=source_type;     // local (0)  gaussian (1)
  header.source_type=wave_state; // S-wave or P-wave source
  header.source_direction=source_direction;
  header.source_laplace_power=LaplacePower;
  header.source_disp_length=disp_length;
  header.source_disp_dir=disp_dir;
  header.sink_smearingparam=0;	// Always to local sinks
  header.sink_type=0;
  header.sink_direction=0;
  header.sink_laplace_power=0;
  header.sink_disp_length=0;
  header.sink_disp_dir=0;

  XMLBufferWriter source_xml;
  write(source_xml, "source", header);

  int ncg_had = 0;

  switch (source_type)
  {
    //
    // Gauge inv. point or shell sources within some S_WAVE, P_WAVE, state etc.
    //
  case SRC_TYPE_POINT_SOURCE:
  case SRC_TYPE_SHELL_SOURCE:
  {
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      
      QDPIO::cout << "color = " << color_source << endl; 

      LatticeColorVector src_color_vec = zero;

      // Make a point source at coordinates t_source
      srcfil(src_color_vec, t_source, color_source);

      // Smear the colour source if specified
      if(source_type == SRC_TYPE_SHELL_SOURCE)
      {

        gausSmear(u_smr, src_color_vec, wvf_param, WvfIntPar, j_decay);
        laplacian(u_smr, src_color_vec, j_decay, LaplacePower);
	displacement(u_smr,src_color_vec,disp_length,disp_dir);

        //power = 1 for one laplacian operator
      }

      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
        QDPIO::cout << "spin = " << spin_source << endl; 

        // Insert a ColorVector into spin index spin_source
        // This only overwrites sections, so need to initialize first
        LatticeFermion chi = zero;

        CvToFerm(src_color_vec, chi, spin_source);
      
        if(wave_state == WAVE_TYPE_P_WAVE)
	  p_src(u_smr, chi, source_direction);

        if(wave_state == WAVE_TYPE_D_WAVE)   /* added */
	  d_src(u_smr, chi, source_direction);

        // primitive initial guess for the linear sys solution
      
        /*
         *  Move the solution to the appropriate components
         *  of quark propagator.
         */

        FermToProp(chi, quark_source, color_source, spin_source);
      }
    }
  }
  break;

  case SRC_TYPE_WALL_SOURCE:
    {
#if 0
    // This stuff is not needed here - should be in propagator.cc

    // For a wall source, we must Coulomb gauge fix the gauge field
    Real GFAccu = 1.0e-6;       // Gauge-fixing relaxation accuracy
    Real OrPara = 1.0;          // Over-relaxation parameter in gauge-fixing
    int GFMax = 1000;           // Maximum number of gauge-fixing relaxations
    bool ORlxDo = false;        // Do Over-relaxation in gauge fixing
    int nrl_gf;                 // Number of relaxations in gauge-fixing

    coulGauge(u, j_decay, GFAccu, GFMax, nrl_gf, ORlxDo, OrPara);
#endif

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
        // Wall fill a fermion source. Insert it into the propagator source
        LatticeFermion chi;
        walfil(chi, t_source[j_decay], j_decay, color_source, spin_source);
        FermToProp(chi, quark_source, color_source, spin_source);
      }
    }
  }
  break;
  
  default:
    QDPIO::cout << "Unsupported source type" << endl;
    QDP_abort(1);
  }
  

  // Sanity check - write out the norm2 of the source in the Nd-1 direction
  // Use this for any possible verification
  {
  // Initialize the slow Fourier transform phases
  SftMom phases(0, true, Nd-1);

  multi1d<Double> source_corr = sumMulti(localNorm2(quark_source),
                                         phases.getSet());

  push(xml_out, "Source_correlator");
  Write(xml_out, source_corr);
  pop(xml_out);
 }
 
  // Now write the source
  switch (prop_type) 
  {
  case PROP_TYPE_SZIN:
    writeSzinQprop(quark_source, source_filename, kappa_fake);
    break;

  case PROP_TYPE_SCIDAC:
    writeQprop(source_filename, quark_source, header);
    break;

  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


#if 0
  // SciDAC output format - move up into switch statement or a subroutine
  {
    XMLBufferWriter file_xml;
    push(file_xml,"make_source");
    write(file_xml, "input", xml_in);
    write(file_xml, "config_info", gauge_xml);
    pop(file_xml);

    QDPSerialFileWriter to(file_xml,source_filename);

    XMLBufferWriter record_xml;
    write(record_xml, "prop_header", header);

    to.write(record_xml,quark_source);  // can keep repeating writes for more records
    to.close();
  }
#endif


  pop(xml_out);  // make_source
  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
