// $Id: t_mesons_w.cc,v 3.1 2006-11-07 21:50:10 edwards Exp $
/*! \file
 *  \brief Test the Wilson mesons() routine
 */

#include <iostream>
#include <cstdio>
#include <time.h>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_mesons_w.xml");
  push(xml,"t_mesons_w");

  push(xml,"lattice");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  pop(xml);

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

  mesons(quark_prop_1, quark_prop_2, phases, t0, xml,
         "Point_Point_Wilson_Mesons") ;

  clock_t end_clock = clock() ;

  QDPIO::cout << "Test took " << end_clock - start_clock << " clocks.\n" ;

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}

