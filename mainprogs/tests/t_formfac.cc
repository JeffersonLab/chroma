// $Id: t_formfac.cc,v 1.1 2003-03-03 16:05:34 flemingg Exp $
/*! \file
 *  \brief Test the form-factor routine
 */

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

  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay];
  multi1d<int> t_source(Nd);
  t_source = 0;

  int t_sink = length-1;

  multi1d<int> sink_mom(Nd-1);

  sink_mom = 0;

  FormFac(u, quark_prop_1, quark_prop_2, t_source, 3, t_sink, sink_mom,
          j_decay, nml);

  // Time to bolt
  QDP_finalize();

  return 0;
}

