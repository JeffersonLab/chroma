// $Id: t_formfac.cc,v 1.3 2003-03-20 19:34:25 flemingg Exp $
//
//! \file
//  \brief Test the form-factor routine
//
// $Log: t_formfac.cc,v $
// Revision 1.3  2003-03-20 19:34:25  flemingg
// Evolved formfac_w.cc to use SftMom class, which included some bug fixes
// in features in SftMom which had been previously untested and evolution
// of the corresponding test program.
//
// Revision 1.2  2003/03/06 00:27:29  flemingg
// Added $Log: t_formfac.cc,v $
// Added Revision 1.3  2003-03-20 19:34:25  flemingg
// Added Evolved formfac_w.cc to use SftMom class, which included some bug fixes
// Added in features in SftMom which had been previously untested and evolution
// Added of the corresponding test program.
// Added to header comments
//

#include <iostream>
#include <cstdio>

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

  NmlWriter nml("t_formfac.nml");

  push(nml,"lattice");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);
  pop(nml);

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

  FormFac(u, quark_prop_1, quark_prop_2, phases, t0, nml) ;

  clock_t end_clock = clock() ;

  cerr << "Test took " << end_clock - start_clock << " clocks.\n" ;

  // Time to bolt
  QDP_finalize();

  return 0;
}

