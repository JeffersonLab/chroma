// $Id: t_dslashm.cc,v 1.1 2002-12-16 07:08:33 edwards Exp $

#include <iostream>
#include <cstdio>

#include <szin.h>

using namespace QDP;


void main_start(void)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the geometry
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::create(nrow);

  //! Test out propagators
  multi1d<LatticeGauge> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;
/* dslash_2d_plus(chi, u, psi, 0); */

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
