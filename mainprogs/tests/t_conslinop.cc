// $Id: t_conslinop.cc,v 1.20 2005-03-02 00:44:19 edwards Exp $

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_conslinop.xml");
  push(xml,"t_conslinop");

  push(xml,"lattis");
  write(xml,"Nd",  Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);

  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Create a fermion BC. Note, the handle is on an ABSTRACT type
  Handle<FermBC<LatticeFermion> >  fbc(new PeriodicFermBC<LatticeFermion>);

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
  
  UnprecWilsonFermAct S(fbc,Mass);
  Handle<const ConnectState> state(S.createState(u));
  Handle<const LinearOperator<LatticeFermion> > A(S.linOp(state));

  LatticeFermion   tmp;
  D(tmp, psi, PLUS);
  DComplex np = innerProduct(psi,tmp);
  D(tmp, psi, MINUS);
  DComplex nm = innerProduct(psi,tmp);

  push(xml,"norm_check");
  write(xml,"np", np);
  write(xml,"nm", nm);
  pop(xml);

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
