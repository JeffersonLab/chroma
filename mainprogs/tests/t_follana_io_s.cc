// $Id: t_follana_io_s.cc,v 1.1 2003-09-11 11:12:21 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "follana_io.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {16,16,16,32};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  // Try and read the propagator;

  LatticePropagator qprop;
  readQpropFollana("./prop.0", qprop);

  NmlWriter nml("prop.nml");

  push(nml, "propagator");
  Write(nml, qprop);
  pop(nml);

  // Time to bolt
  QDP_finalize();
}
