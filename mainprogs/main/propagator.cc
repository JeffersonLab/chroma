// $Id: propagator.cc,v 1.15 2003-06-19 17:34:25 ikuro Exp $
/*! \file
 *  \brief Main code for propagator generation
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
#define D_WAVE 2    /*added*/

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
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  // Useful parameters that should be read from an input file
  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay]; // define the temporal direction


  /*
   *  As a temporary measure, we will now read in the parameters from a file
   *  DATA.  Eventually, this will use QIO or NML, but for the moment we will
   *  just use the usual command-line reader
   */

  TextReader params_in("DATA");


  int io_version_in; 		// The I/O version that we are reading....

  Real Kappa;			// Kappa value
  
  int source_type, source_direction; // S-wave(0), P-wave(1), D-wave(2), and direction

  int wf_type;			// Point (0) or Smeared (1)
  Real wvf_param;		// Parameter for the wave function
  int WvfIntPar;

  Real RsdCG;
  int MaxCG;			// Iteration parameters

  params_in >> io_version_in;

  switch(io_version_in){	// The parameters we read in IO version

  case 101:			// 

    params_in >> Kappa;

    params_in >> source_type;	// S-wave, P-wave D-wave, etc
    params_in >> source_direction; // dx(0) dy(1) dz(2) dydz(3) dzdz(4)

    params_in >> wf_type;	// Point, Gaussian(2) etc
    params_in >> wvf_param;
    params_in >> WvfIntPar;

    params_in >> RsdCG;		// Target residue and maximum iterations
    params_in >> MaxCG;
    break;

  default:

    QDP_error_exit("Unknown io version", io_version_in);

  }
  
  cout << "Kappa is " << Kappa << endl;

  switch(wf_type){
  case POINT_SOURCE:  // 0
    cout << "Point source" << endl;
    break;
  case SHELL_SOURCE:  // 2
    cout << "Smeared source wvf_param= " << wvf_param <<": WvfIntPar= " 
	 << WvfIntPar << endl;
    break;
  default:
    QDP_error_exit("Unknown source_type", wf_type);
  }    

  cout << "RsdCG= " << RsdCG << ": MaxCG= " << MaxCG << endl;

  multi1d<int> t_source(Nd);
  t_source = 0;


  UnprecWilsonFermAct S_f(Kappa);

  FermAct = UNPRECONDITIONED_WILSON;  // global
  InvType = CG_INVERTER;  // global

  cerr << "DEBUG 1" << endl;

  // Generate a hot start gauge field
  multi1d<LatticeColorMatrix> u(Nd);

  cerr << "DEBUG 1" << endl;


  /*  for(int mu=0; mu < u.size(); ++mu)
  {
    gaussian(u[mu]);
    reunit(u[mu]);
    }
  */

  readArchiv(u, "nersc_freefield.cfg");

  //Seed seed_old;
  //readSzin(u, "szin.cfg", seed_old);

  cerr << "DEBUG 2" << endl;


  // Useful info

  string nml_filename;

  switch(source_type){
  case S_WAVE:
    nml_filename = "propagator.nml";
    break;
  case P_WAVE:
    nml_filename = "dz_propagator.nml";
    break;
  case D_WAVE:    /* added */
    if (source_direction == 12)
      nml_filename = "dydz_propagator.nml";
    if (source_direction == 22)
      nml_filename = "dzdz_propagator.nml";
    break;
  default: 
    cerr<<"invaid source_type\n";
    break;
  }

  cerr << "DEBUG 3" << endl;

  NmlWriter nml(nml_filename);
    

  push(nml,"lattice");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);
  pop(nml);

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
  header.kappa = Kappa;
  header.source_smearingparam=wf_type;     // local (0)  gaussian (2)
  header.source_type=source_type; // S-wave, P-wave or D-wave source
  header.source_direction=source_direction; 
  header.sink_smearingparam=0;	// Always to local sinks
  header.sink_type=0;
  header.sink_direction=0;   // dx(0) dy(1) dz(2) dydy(11) dydz(12) dzdz(22)
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

    if(wf_type == SHELL_SOURCE) {
      gausSmear(u, src_color_vec, wvf_param, WvfIntPar, j_decay);
      //laplacian(u, src_color_vec, j_decay, power); 
      // power = 1 for one laplacian operator
      //         2 for two 
    }


    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      if (Layout::primaryNode())
        cerr << "spin = " << spin_source << endl;

      // Insert a ColorVector into spin index spin_source
      // This only overwrites sections, so need to initialize first

      LatticeFermion psi = zero;  // note this is ``zero'' and not 0

      {
	
	LatticeFermion chi = zero;

	CvToFerm(src_color_vec, chi, spin_source);

	if(source_type == P_WAVE)
	  p_src(u, chi, source_direction);

	if(source_type == D_WAVE)   /* added */
	  d_src(u, chi, source_direction);



	// primitive initial guess for the linear sys solution

	// Compute the propagator for given source color/spin.
	push(nml,"qprop");
	int n_count;

	S_f.qprop(psi, u, chi, RsdCG, MaxCG, n_count);
	ncg_had += n_count;
	
	Write(nml, Kappa);
	Write(nml, RsdCG);
	Write(nml, n_count);
	  
	pop(nml);

	/*
	 *  Move the solution to the appropriate components
	 *  of quark propagator.
	 */
      }
      FermToProp(psi, quark_propagator, color_source, spin_source);
    }
  }

  switch(source_type){
  case S_WAVE:
    writeQprop("propagator_0", quark_propagator, header);
    break;
  case P_WAVE:
    writeQprop("dz_propagator_0", quark_propagator, header);
    break;
  case D_WAVE:       /* added */
    if (source_direction ==12)
      writeQprop("dydz_propagator_0", quark_propagator, header);
    if (source_direction ==22)
      writeQprop("dzdz_propagator_0", quark_propagator, header);
    break;
  default:
    QDP_error_exit("Unknown io version", io_version_in);
  }    

  //Write(nml, quark_propagator);

  nml.close();

  // Time to bolt
  QDP_finalize();

  return 0;
}
