// $Id: t_dwflinop.cc,v 1.1 2003-11-08 04:18:59 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

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

  XMLFileWriter xml("t_dwflinop.xml");
  push(xml, "t_dwflinop");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  Real WilsonMass = 1.5;
  Real m_q = 0.1;
  UnprecDWFermAct S_f(WilsonMass, m_q);

  const  LinearOperator<LatticeDWFermion>* A = S_f.linOp(u);

  LatticeDWFermion psi, chi;
  random(psi);
  random(chi);

  chi = (*A)(psi, PLUS);
  DComplex nn1 = innerProduct(psi, chi);

  chi = (*A)(psi, MINUS);
  DComplex nn2 = innerProduct(chi, psi);

  push(xml,"innerprods");
  Write(xml, nn1);
  Write(xml, nn2);
  pop(xml);

  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
