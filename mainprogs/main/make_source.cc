// $Id: make_source.cc,v 1.1 2003-06-20 20:47:35 dgr Exp $
/*! \file
 *  \brief Main code for source generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

/*
 *  Here we have various temporary definitions
 */

// First the source type
#define S_WAVE 0
#define P_WAVE 1

#define MAXLINE 80

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

  // Setup the layout
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);

  // Useful parameters that should be read from an input file
  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay]; // define the temporal direction


  /*
   *  As a temporary measure, we will now read in the parameters from a file
   *  DATA.  Eventually, this will use QIO or NML, but for the moment we will
   *  just use the usual command-line reader
   */

  NmlReader nml_in("DATA");


  int version; 		// The input-parameter version

  Real kappa_fake = 0.0;			// Kappa value
  
  int source_type, source_direction; // S-wave, P-wave etc, and direction

  int wf_type;			// Point (0) or Smeared (1)
  Real wvf_param;		// Parameter for the wave function
  int WvfIntPar;

  Real RsdCG;
  int MaxCG;			// Iteration parameters

  push(nml_in, "IO_version") ;
  Read(nml_in, version) ;
  pop(nml_in) ;


  switch(version){	// The parameters we read in IO version

  case 101:			

    push(nml_in,"param");	// Push into param group

    Read(nml_in, source_type);	// S-wave, P-wave etc
    Read(nml_in, source_direction);

    Read(nml_in, wf_type);	// Point, Gaussian etc
    Read(nml_in,  wvf_param);
    Read(nml_in, WvfIntPar);

    break;

  default:

    QDP_error_exit("Unknown io version", version);

  }
  
  // Now get the lattice sizes etc
  Read(nml_in, nrow);


  nml_in.close();

  switch(wf_type){
  case POINT_SOURCE:
    cout << "Point source" << endl;
    break;
  case SHELL_SOURCE:
    cout << "Smeared source wvf_param= " << wvf_param <<": WvfIntPar= " 
	 << WvfIntPar << endl;
    break;
  default:
    QDP_error_exit("Unknown source_type", wf_type);
  }    

  Layout::setLattSize(nrow);
  Layout::create();



  multi1d<int> t_source(Nd);
  t_source = 0;

  // Generate a hot start gauge field
  multi1d<LatticeColorMatrix> u(Nd);


  /*  readArchiv(u, "nersc_freefield.cfg");*/

  Seed seed_old;
  readSzin(u, "szin.cfg", seed_old);


  // Useful info

  string nml_filename;
  string source_filename;


  switch(source_type){
  case S_WAVE:
    nml_filename = "source.nml";
    source_filename = "source_0";
    break;
  case P_WAVE:
    nml_filename = "p_source.nml";
    source_filename = "p_source_0";

  }


  NmlWriter nml_out(nml_filename);
    

  push(nml_out,"lattice");
  Write(nml_out,Nd);
  Write(nml_out,Nc);
  Write(nml_out,Ns);
  Write(nml_out,nrow);
  pop(nml_out);

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

#if 1
  PropHead header;		// Header information
  header.kappa = kappa_fake;
  header.source_smearingparam=wf_type;     // local (0)  gaussian (1)
  header.source_type=source_type; // S-wave or P-wave source
  header.source_direction=source_direction;
  header.sink_smearingparam=0;	// Always to local sinks
  header.sink_type=0;
  header.sink_direction=0;
#endif

  int ncg_had = 0;


  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    if (Layout::primaryNode())
      cerr << "color = " << color_source << endl;

    LatticeColorVector src_color_vec = zero;

    // Make a point source at coordinates t_source

    srcfil(src_color_vec, t_source, color_source);

    // Smear the colour source if specified

    if(wf_type == SHELL_SOURCE)
      gausSmear(u, src_color_vec, wvf_param, WvfIntPar, j_decay);


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


      // primitive initial guess for the linear sys solution
      
      /*
       *  Move the solution to the appropriate components
       *  of quark propagator.
       */

      FermToProp(chi, quark_propagator, color_source, spin_source);
    }
  }


  writeSzinQprop(quark_propagator, source_filename, kappa_fake);


  nml_out.close();

  // Time to bolt
  QDP_finalize();

  return 0;
}
