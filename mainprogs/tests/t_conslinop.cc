// $Id: t_conslinop.cc,v 1.3 2003-04-03 19:28:07 edwards Exp $

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

  NmlWriter nml("t_conslinop.nml");

  push(nml,"lattis");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,nrow);
  pop(nml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);

  cerr << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);


  Real Kappa = 0.1;
//  LinearOperator* A = ConsLinOp(u, Kappa, UNPRECONDITIONED_WILSON);
  WilsonDslash D(u);
//  WilsonDslash D;

  LatticeFermion psi, chi;
  gaussian(psi);
  cerr << "before dslash call" << endl;
  chi[rb[0]] = D(psi, PLUS, 0); 
  chi[rb[1]] = D(psi, PLUS, 0); 
  cerr << "after dslash call" << endl;

  cerr << "before wilson construct" << endl;
  UnpreconditionedWilson M(u,Kappa);
  cerr << "after wilson construct" << endl;
  chi[rb[1]] = M(psi, PLUS); 
  cerr << "after wilson call" << endl;
  
  // Time to bolt
  QDP_finalize();
}
