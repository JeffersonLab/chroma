// $Id: t_conslinop.cc,v 1.13 2003-11-20 05:43:41 edwards Exp $

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

  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);


  Real Mass = 0.1;
  WilsonDslash D(u);
//  WilsonDslash D;

  LatticeFermion psi, chi;
  gaussian(psi);
  QDPIO::cout << "before dslash call" << endl;
  D.apply(chi, psi, PLUS, 0); 
  D.apply(chi, psi, PLUS, 1); 
  QDPIO::cout << "after dslash call" << endl;

  QDPIO::cout << "before wilson construct" << endl;
  UnprecWilsonLinOp M(u,Mass);
  QDPIO::cout << "after wilson construct" << endl;
  M(chi, psi, PLUS); 
  QDPIO::cout << "after wilson call" << endl;
  
  UnprecWilsonFermAct S(Mass);
  const LinearOperator<LatticeFermion>* A = S.linOp(u);

  LatticeFermion   tmp;
  D(tmp, psi, PLUS);
  DComplex np = innerProduct(psi,tmp);
  D(tmp, psi, MINUS);
  DComplex nm = innerProduct(psi,tmp);

  push(nml,"norm_check");
  Write(nml,np);
  Write(nml,nm);
  pop(nml);

  delete A;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
