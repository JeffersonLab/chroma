/*
 *  $Id: qqq_w.cc,v 1.3 2003-06-19 17:34:25 ikuro Exp $
 *
 *  This is the test program for the routine that reads in a quark propagator,
 *  stored in SZIN format, and computes the generalised quark propagators
 *  that will be the basis of our baryon code
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace QDP;

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

  multi1d<LatticeColorMatrix> u(Nd);
  // readArchiv(u, "nersc_freefield.cfg");	

  Seed seed_old;
  readSzin(u, "szin.cfg", seed_old);
  


  // Now the lattice quark propagator, just a single one for this example, together
  // with the corresponding header

  LatticePropagator s_s_quark_propagator;
  PropHead          s_s_header;
  readQprop("propagator_0", s_s_quark_propagator, s_s_header);

  NmlWriter nml("GP.nml");

  push(nml, "qqqqqq_barcomp");
  barcomp(s_s_quark_propagator, s_s_header,
	  s_s_quark_propagator, s_s_header,
	  s_s_quark_propagator, s_s_header,
	  0, j_decay, 1, "qqqqqq_barcomp", nml);
      // t0, j_decay, bc


  /*******************************************************************/
  /***************************** P-WAVE ******************************/
  /*******************************************************************/

  {
    LatticePropagator s_p_quark_propagator;
    PropHead          s_p_header;
    readQprop("dz_propagator_0", s_p_quark_propagator, s_p_header); 

    {
      LatticePropagator p_p_quark_propagator;
      PropHead          p_p_header = s_p_header;

      push(nml, "qqqqDzqDzq_barcomp");
      p_p_header.sink_type=1;
      p_p_header.sink_direction=2;  
      D_j(u,  s_p_quark_propagator, p_p_quark_propagator, 2);
      barcomp(s_s_quark_propagator, s_s_header,
	      s_s_quark_propagator, s_s_header,
	      p_p_quark_propagator, p_p_header,
	      0, j_decay, 1, "qqqqDzqDzq_barcomp", nml);
      
      push(nml, "qqqqDxqDzq_barcomp");
      p_p_header.sink_type=1;
      p_p_header.sink_direction=0;  
      D_j(u,  s_p_quark_propagator, p_p_quark_propagator, 0);
      barcomp(s_s_quark_propagator, s_s_header,
	      s_s_quark_propagator, s_s_header,
	      p_p_quark_propagator, p_p_header,
	      0, j_decay, 1, "qqqqDxqDzq_barcomp", nml);
      
      push(nml, "qqqqDyqDzq_barcomp");
      p_p_header.sink_type=1;
      p_p_header.sink_direction=1;  
      D_j(u,  s_p_quark_propagator, p_p_quark_propagator, 1);
      barcomp(s_s_quark_propagator, s_s_header,
	      s_s_quark_propagator, s_s_header,
	      p_p_quark_propagator, p_p_header,
	      0, j_decay, 1, "qqqqDyqDzq_barcomp", nml);
    }

    {  
      LatticePropagator p_s_quark_propagator;
      PropHead          p_s_header = s_s_header;

      push(nml, "qqDzqqqDzq_barcomp");
      p_s_header.sink_type=1;
      p_s_header.sink_direction=2;  
      D_j(u,  s_s_quark_propagator, p_s_quark_propagator, 2);
      barcomp(s_s_quark_propagator, s_s_header,
	      p_s_quark_propagator, p_s_header,
	      s_p_quark_propagator, s_p_header,
	      0, j_decay, 1, "qqDzqqqDzq_barcomp", nml);
      
      push(nml, "qqDxqqqDzq_barcomp");
      p_s_header.sink_type=1;
      p_s_header.sink_direction=0;  
      D_j(u,  s_s_quark_propagator, p_s_quark_propagator, 0);
      barcomp(s_s_quark_propagator, s_s_header,
	      p_s_quark_propagator, p_s_header,
	      s_p_quark_propagator, s_p_header,
	      0, j_decay, 1, "qqDxqqqDzq_barcomp", nml);
      
      push(nml, "qqDyqqqDzq_barcomp");
      p_s_header.sink_type=1;
      p_s_header.sink_direction=1;  
      D_j(u,  s_s_quark_propagator, p_s_quark_propagator, 1);
      barcomp(s_s_quark_propagator, s_s_header,
	      p_s_quark_propagator, p_s_header,
	      s_p_quark_propagator, s_p_header,
	      0, j_decay, 1, "qqDyqqqDzq_barcomp", nml);
    }
  }

  /*******************************************************************/
  /***************************** D-WAVE ******************************/
  /*******************************************************************/

  {
    LatticePropagator s_dydz_quark_propagator;
    PropHead          s_dydz_header;
    readQprop("dydz_propagator_0", s_dydz_quark_propagator, s_dydz_header);
    // dydz_propagator_0 contains s_dydz_quark_propagators

    {
      LatticePropagator dydz_dydz_quark_propagator;
      PropHead          dydz_dydz_header = s_dydz_header;

      push(nml, "qqqqDyDzqDyDzq_barcomp");
      dydz_dydz_header.sink_type=2;
      dydz_dydz_header.sink_direction=12; // mean dydz
      DjDk(u, s_dydz_quark_propagator, dydz_dydz_quark_propagator, 12);
      //                                                    1->dy 2->dz
      barcomp(s_s_quark_propagator,       s_s_header,
	      s_s_quark_propagator,       s_s_header,
	      dydz_dydz_quark_propagator, dydz_dydz_header,
	      0, j_decay, 1, "qqqqDyDzqDyDzq_barcomp", nml);
    }

    {
      LatticePropagator dydz_s_quark_propagator;
      PropHead          dydz_s_header = s_s_header;

      push(nml, "qqDyDzqqqDyDzq_barcomp");
      dydz_s_header.sink_type=2;
      dydz_s_header.sink_direction=12;  // mean dydz
      DjDk(u, s_s_quark_propagator, dydz_s_quark_propagator, 12);
      barcomp(s_s_quark_propagator,    s_s_header,
	      dydz_s_quark_propagator, dydz_s_header,
	      s_dydz_quark_propagator, s_dydz_header,
	      0, j_decay, 1, "qqDyDzqqqDyDzq_barcomp", nml);
    }
  }


  {
    LatticePropagator s_dzdz_quark_propagator;
    PropHead          s_dzdz_header;
    readQprop("dzdz_propagator_0", s_dzdz_quark_propagator, s_dzdz_header);
    // dzdz_propagator_0 contains s_dzdz_quark_propagators

    {
      LatticePropagator dzdz_dzdz_quark_propagator;
      PropHead          dzdz_dzdz_header = s_dzdz_header;

      push(nml, "qqqqDzDzqDzDzq_barcomp");
      dzdz_dzdz_header.sink_type=2;
      dzdz_dzdz_header.sink_direction=22; // mean dzdz
      DjDk(u, s_dzdz_quark_propagator, dzdz_dzdz_quark_propagator, 22);
      barcomp(s_s_quark_propagator,       s_s_header,
	      s_s_quark_propagator,       s_s_header,
	      dzdz_dzdz_quark_propagator, dzdz_dzdz_header,
	      0, j_decay, 1, "qqqqDzDzqDzDzq_barcomp", nml);
    }

    {
      LatticePropagator dzdz_s_quark_propagator;
      PropHead          dzdz_s_header = s_s_header;

      push(nml, "qqDzDzqqqDzDzq_barcomp");
      dzdz_s_header.sink_type=2;
      dzdz_s_header.sink_direction=22;
      DjDk(u, s_s_quark_propagator, dzdz_s_quark_propagator, 22);
      barcomp(s_s_quark_propagator,    s_s_header,
	      dzdz_s_quark_propagator, dzdz_s_header,
	      s_dzdz_quark_propagator, s_dzdz_header,
	      0, j_decay, 1, "qqDzDzqqqDzDzq_barcomp", nml);
    }

    {
      LatticePropagator dydy_dzdz_quark_propagator;
      PropHead          dydy_dzdz_header = s_dzdz_header;

      push(nml, "qqqqDyDyqDzDzq_barcomp");
      dydy_dzdz_header.sink_type=2;
      dydy_dzdz_header.sink_direction=11;  // mean dydy
      DjDk(u, s_dzdz_quark_propagator, dydy_dzdz_quark_propagator, 11);
      barcomp(s_s_quark_propagator,       s_s_header,
	      s_s_quark_propagator,       s_s_header,
	      dydy_dzdz_quark_propagator, dydy_dzdz_header,
	      0, j_decay, 1, "qqqqDyDyqDzDzq_barcomp", nml);
    }
      
    {
      LatticePropagator dydy_s_quark_propagator;
      PropHead          dydy_s_header = s_s_header;

      push(nml, "qqDyDyqqqDzDzq_barcomp");
      dydy_s_header.sink_type=2;
      dydy_s_header.sink_direction=11;  // mean dydy
      DjDk(u, s_s_quark_propagator, dydy_s_quark_propagator, 11);
      barcomp(s_s_quark_propagator,    s_s_header,
	      dydy_s_quark_propagator, dydy_s_header,
	      s_dzdz_quark_propagator, s_dzdz_header,
	      0, j_decay, 1, "qqDyDyqqqDzDzq_barcomp", nml);
    }
  }




  nml.close();
}
