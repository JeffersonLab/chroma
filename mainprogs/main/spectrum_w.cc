// $Id: spectrum_w.cc,v 1.1 2003-04-10 20:20:18 dgr Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace QDP;

/*
 *  This routine computes the meson and baryon spectrum for
 *  a single propagator, propagator_0
 *
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

  multi1d<int> t_source(Nd);
  t_source = 0;

  SftMom phases(0, true, j_decay); // Set up the momentum phases

  LatticePropagator quark_prop;
  PropHead quark_head;

  NmlWriter nml("spectrum.nml");

  push(nml,"lattice");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);
  pop(nml);

  readQprop("propagator_0", quark_prop, quark_head);

  mesons(quark_prop, quark_prop, phases, t_source[j_decay], nml,
         "Point_Point_Wilson_Mesons") ;


  nml.close();

  // Time to bolt
  QDP_finalize();

  return 0;
}
