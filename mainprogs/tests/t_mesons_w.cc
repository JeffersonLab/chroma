// $Id: t_mesons_w.cc,v 1.2 2003-03-14 05:14:32 flemingg Exp $
//
//! \file
//  \brief Test the Wilson mesons() routine
//
// $Log: t_mesons_w.cc,v $
// Revision 1.2  2003-03-14 05:14:32  flemingg
// rewrite of mesons_w.cc to use the new SftMom class.  mesons_w.cc still
// needs to be cleaned up once the best strategy is resolved.  But for now,
// the library and test program compiles and runs.
//

#include <iostream>
#include <cstdio>
#include <time.h>

#include "chroma.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  NmlWriter nml("t_mesons_w.nml");

  push(nml,"lattice");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);
  pop(nml);

  LatticePropagator quark_prop_1, quark_prop_2;
  gaussian(quark_prop_1);
  gaussian(quark_prop_2);

  // propagation direction
  int j_decay = Nd-1 ;

  // source timeslice
  int t0 = 0 ;

  clock_t start_clock = clock() ;
  // create averaged Fourier phases with (mom)^2 <= 10
  SftMom phases(10, true, j_decay) ;

  mesons(quark_prop_1, quark_prop_2, phases, t0, nml,
         "Point_Point_Wilson_Mesons") ;

  clock_t end_clock = clock() ;

  cerr << "Test took " << end_clock - start_clock << " clocks.\n" ;

  // Time to bolt
  QDP_finalize();

  return 0;
}

