// $Id: t_conslinop.cc,v 1.11 2003-09-11 00:46:04 edwards Exp $

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
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
  InvType = CG_INVERTER;
  WilsonDslash D(u);
//  WilsonDslash D;

  LatticeFermion psi, chi;
  gaussian(psi);
  cerr << "before dslash call" << endl;
  chi[rb[0]] = D.apply(psi, PLUS, 0); 
  chi[rb[1]] = D.apply(psi, PLUS, 0); 
  cerr << "after dslash call" << endl;

  cerr << "before wilson construct" << endl;
  UnprecWilsonLinOp M(u,Kappa);
  cerr << "after wilson construct" << endl;
  chi = M(psi, PLUS); 
  cerr << "after wilson call" << endl;
  
  UnprecWilsonFermAct S(Kappa);
  const LinearOperator* A = S.linOp(u);

  DComplex np = innerProduct(psi,D(psi,PLUS));
  DComplex nm = innerProduct(psi,D(psi,MINUS));

  push(nml,"norm_check");
  Write(nml,np);
  Write(nml,nm);
  pop(nml);

  delete A;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
