// $Id: t_dslashm.cc,v 1.9 2004-11-22 22:54:39 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;
  dslash(chi, u, psi, +1, 0);

  XMLFileWriter xml("t_dslashm.xml");
  push(xml,"t_dslashm");

  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  write(xml,"psi", psi);
  write(xml,"chi", chi);

  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
