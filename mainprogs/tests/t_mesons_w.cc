// $Id: t_mesons_w.cc,v 1.1 2003-03-06 00:30:15 flemingg Exp $
//
//! \file
//  \brief Test the Wilson mesons() routine
//
// $Log: t_mesons_w.cc,v $
// Revision 1.1  2003-03-06 00:30:15  flemingg
// Complete rewrite of lib/meas/hadron/mesons_w.cc, including a simple test
// program in mainprogs/tests built with 'make check' and various other
// changes to autoconf/make files to support this rewrite.
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

  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay];
  multi1d<int> t_source(Nd);
  t_source = 0;

  mesons(quark_prop_1, quark_prop_2, t_source, 10, j_decay, nml);

  // Time to bolt
  QDP_finalize();

  return 0;
}

