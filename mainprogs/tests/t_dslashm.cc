// $Id: t_dslashm.cc,v 1.2 2003-01-04 05:09:27 edwards Exp $

#include <iostream>
#include <cstdio>

#include "szin.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::create(nrow);

  //! Test out propagators
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;
  dslash(chi, u, psi, +1, 0);

  NmlWriter nml("t_dslashm.nml");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);
  Write(nml,psi);
  Write(nml,chi);

  // Time to bolt
  QDP_finalize();
}
