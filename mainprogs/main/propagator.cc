// $Id: propagator.cc,v 1.5 2003-04-04 18:09:06 edwards Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace QDP;

//! Main propagator generation code
/*! \defgroup propagator Main code for propagator generation
 *  \ingroup main
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  // Useful parameters that should be read from an input file
  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay]; // define the temporal direction

  multi1d<int> t_source(Nd);
  t_source = 0;

  Real wvf_param = 1.1;  // smearing width
  int  WvfIntPar = 10;   // number of hits for smearing

  Real Kappa = 0.1480;

  FermAct = UNPRECONDITIONED_WILSON;  // global
  InvType = CG_INVERTER;  // global
  Real RsdCG = 1.0e-6;
  int  MaxCG = 500;

  // Generate a hot start gauge field
  multi1d<LatticeColorMatrix> u(Nd);

  for(int mu=0; mu < u.size(); ++mu)
  {
    gaussian(u[mu]);
    reunit(u[mu]);
  }

  /*  readArchiv(u, "nersc_freefield.cfg");*/

  //  Seed seed_old;
  //  readSzin(u, 0, "szin.cfg", seed_old);

  // Useful info
  NmlWriter nml("propagator.nml");

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

#if 0
  PropHead header;		// Header information
  header.kappa = Kappa;
  header.source_smearingparam=0;     // local (0)  gaussian (1)
  header.source_type=0;		// S-wave source
  header.source_direction=0;
  header.sink_smearingparam=0;
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

    // Smear the color source
    //    gausSmear(u, src_color_vec, wvf_param, WvfIntPar, j_decay);	

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

	cerr << "DEBUG After CvToFerm " << endl;


	// primitive initial guess for the linear sys solution

	// Compute the propagator for given source color/spin.
	push(nml,"Qprop");
	int n_count;

	cerr << "DEBUG before NPQprop " << endl;
	NPQprop(u, chi, Kappa, RsdCG, psi, MaxCG, n_count);
	cerr << "DEBUG after NPQprop " << endl;
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

//  writeQprop("propagator_0", quark_propagator, header);
  BinaryWriter cfg_out("propagator_0");
  write(cfg_out,quark_propagator);
  cfg_out.close();

  nml.close();

  // Time to bolt
  QDP_finalize();

  return 0;
}
