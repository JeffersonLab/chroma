// $Id: make_source.cc,v 1.8 2003-10-01 02:51:25 edwards Exp $
/*! \file
 *  \brief Main code for source generation
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"

/*
 *  Here we have various temporary definitions
 */
enum CfgType {
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_UNKNOWN
} ;

enum PropType {
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
  PROP_TYPE_UNKNOWN
} ;


// First the source type
#define S_WAVE 0
#define P_WAVE 1
#define D_WAVE 2    /*added*/

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

  Real kappa_fake = 0.0;			// Kappa value
  
  int source_type, source_direction; // S-wave(0), P-wave(1), D-wave(2)and direction

  int wf_type;			// Point (0) or Smeared (2)
  Real wvf_param;		// smearing width
  int WvfIntPar;                // number of iteration for smearing
  int LaplacePower;             // power of laplacian operator: 0, 1, 2, etc.

  multi1d<int> t_source(Nd);

  CfgType cfg_type;
  PropType prop_type;

  string cfg_file;

  string xml_in_root = "/make_source";
  string path = xml_in_root + "/IO_version"; // push into 'IO_version' group

  try {
  read(xml_in, path + "/version", version) ;

  switch(version){	// The parameters we read in IO version

  case 2:			

    path = xml_in_root + "/param"; // push into 'param' group

    read(xml_in, path + "/j_decay", j_decay);

    read(xml_in, path + "/source_type", source_type);	// S-wave, P-wave etc
    read(xml_in, path + "/source_direction", source_direction);

    read(xml_in, path + "/wf_type", wf_type);	// Point, Gaussian etc
    read(xml_in, path + "/wvf_param", wvf_param);
    read(xml_in, path + "/WvfIntPar", WvfIntPar);
    read(xml_in, path + "/LaplacePower", LaplacePower);
    read(xml_in, path + "/t_srce", t_source);

    {
      string input_cfg_type ;
      read(xml_in, path + "/cfg_type", input_cfg_type) ;
      if (input_cfg_type == "SZIN")
        cfg_type = CFG_TYPE_SZIN ;
      else 
        cfg_type = CFG_TYPE_UNKNOWN ;
    }
    {
      string input_prop_type ;
      read(xml_in, path + "/prop_type", input_prop_type) ;
      if (input_prop_type == "SZIN")
        prop_type = PROP_TYPE_SZIN ;
      else if (input_prop_type == "SCIDAC")
        prop_type = PROP_TYPE_SCIDAC ;
      else 
        prop_type = PROP_TYPE_UNKNOWN ;
    }
    break;

  default:

    QDP_error_exit("Unknown io version", version);

  }
  
  // Now get the lattice sizes etc
  read(xml_in, path + "/nrow", nrow);

  // Read in the gauge configuration file name
  read(xml_in, xml_in_root + "/Cfg", cfg_file);

  } 
  catch(const string& e)
  {
    QDP_error_exit("Error reading in make_source: %s", e.c_str());
  }

  Layout::setLattSize(nrow);
  Layout::create();

  switch(wf_type){
  case POINT_SOURCE:
    cout << "Point source" << endl;
    break;
  case SHELL_SOURCE:
    cout << "Smeared source wvf_param= " << wvf_param <<": WvfIntPar= " 
	 << WvfIntPar << endl
         << "Power of Laplacian operator" << LaplacePower << endl;
    break;
  default:
    QDP_error_exit("Unknown source_type", wf_type);
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

  // Set the rng seed
  QDP::RNG::setrn (seed_old);


  // Useful info

  string xml_filename;
  string source_filename;


  switch(source_type){
  case S_WAVE:
    xml_filename = "source.xml";
    source_filename = "source_0";
    break;
  case P_WAVE:
    xml_filename = "p_source.xml";
    source_filename = "p_source_0";
    break;
  case D_WAVE:    /* added */
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
    cerr<<"invaid source_type\n";
    break;
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
  LatticePropagator quark_propagator;

  PropHead header;		// Header information
  header.kappa = kappa_fake;
  header.source_smearingparam=wf_type;     // local (0)  gaussian (1)
  header.source_type=source_type; // S-wave or P-wave source
  header.source_direction=source_direction;
  header.source_laplace_power=LaplacePower;
  header.sink_smearingparam=0;	// Always to local sinks
  header.sink_type=0;
  header.sink_direction=0;
  header.sink_laplace_power=0;

  XMLBufferWriter source_xml;
  write(source_xml, "source", header);

  int ncg_had = 0;


  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    if (Layout::primaryNode())
      cerr << "color = " << color_source << endl;

    LatticeColorVector src_color_vec = zero;

    // Make a point source at coordinates t_source

    srcfil(src_color_vec, t_source, color_source);

    // Smear the colour source if specified

    if(wf_type == SHELL_SOURCE) {
      gausSmear(u, src_color_vec, wvf_param, WvfIntPar, j_decay);
      laplacian(u, src_color_vec, j_decay, LaplacePower);
      //power = 1 for one laplacian operator
    }

    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      if (Layout::primaryNode())
        cerr << "spin = " << spin_source << endl;

      // Insert a ColorVector into spin index spin_source
      // This only overwrites sections, so need to initialize first


      LatticeFermion chi = zero;

      CvToFerm(src_color_vec, chi, spin_source);
      
      if(source_type == P_WAVE)
	p_src(u, chi, source_direction);

      if(source_type == D_WAVE)   /* added */
	d_src(u, chi, source_direction);

      // primitive initial guess for the linear sys solution
      
      /*
       *  Move the solution to the appropriate components
       *  of quark propagator.
       */

      FermToProp(chi, quark_propagator, color_source, spin_source);
    }
  }


  switch (prop_type) 
  {
  case PROP_TYPE_SZIN:
    writeSzinQprop(quark_propagator, source_filename, kappa_fake);
    break;

  case PROP_TYPE_SCIDAC:
    writeQprop(source_filename, quark_propagator, header);
    break;

  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


#if 1
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

    to.write(record_xml,quark_propagator);  // can keep repeating writes for more records
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
