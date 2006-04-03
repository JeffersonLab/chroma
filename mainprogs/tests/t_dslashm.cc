// $Id: t_dslashm.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

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
  Chroma::finalize();

  exit(0);
}
