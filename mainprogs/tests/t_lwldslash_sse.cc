// $Id: t_lwldslash_sse.cc,v 1.1 2003-09-10 18:15:06 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "primitives.h" // GTF: for PLUS
#include "lib.h"


using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  NmlWriter nml("t_lwldslash.nml");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi1, chi2;
  random(psi);
  chi1 = zero;
  chi2 = zero;

  //! Create a linear operator
  WilsonDslash D_w(u);
  /* SSEWilsonDslash D_sse(u); */

  chi1 = D_w.apply(psi, PLUS, 0);
  /* chi2 = D_sse.apply(psi,PLUS, 0); */

  LatticeFermion res= chi2 - chi1;

  Real n = norm2(res);
  cout << "Norm of diff is " << n;
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);

  // Time to bolt
  QDP_finalize();

  return 0;
}
