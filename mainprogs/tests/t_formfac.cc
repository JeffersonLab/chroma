// $Id: t_formfac.cc,v 3.1 2006-11-07 21:50:10 edwards Exp $
/*! \file
 *  \brief Test the form-factor routine
 */

#include <iostream>
#include <cstdio>

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

  XMLFileWriter xml("t_formfac.xml");

  push(xml,"lattice");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  pop(xml);

  // Randomize the gauge field
  multi1d<LatticeColorMatrix> u(Nd);

  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

  LatticePropagator quark_prop_1, quark_prop_2;
  gaussian(quark_prop_1);
  gaussian(quark_prop_2);

  // propagation direction
  int j_decay = Nd-1;

  // source timeslice
  int t0 = 0 ;

  // sink momentum {1,0,0}
  multi1d<int> sink_mom(Nd-1) ;
  sink_mom    = 0 ;
  sink_mom[0] = 1 ;

  clock_t start_clock = clock() ;

  // create Fourier phases with (mom-sink_mom)^2 <= 10 (NO AVERAGING)
  SftMom phases(10, sink_mom, false, j_decay) ;

  FormFac(u, quark_prop_1, quark_prop_2, phases, t0, xml) ;

  clock_t end_clock = clock() ;

  QDPIO::cout << "Test took " << end_clock - start_clock << " clocks.\n" ;

  // Time to bolt
  Chroma::finalize();

  exit(0);
}

