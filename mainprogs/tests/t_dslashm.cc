// $Id: t_dslashm.cc,v 1.8 2004-02-11 12:51:35 bjoo Exp $

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

  NmlWriter nml("t_dslashm.nml");
  write(nml,"Nd", Nd);
  write(nml,"Nc", Nc);
  write(nml,"Ns", Ns);
  write(nml,"nrow", nrow);
  write(nml,"psi", psi);
  write(nml,"chi", chi);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
