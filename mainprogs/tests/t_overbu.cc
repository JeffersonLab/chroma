// $Id: t_overbu.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  //! Test out overbu
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
  {
//    gaussian(u[m]);
    u[m] = 1.0;
  }

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  gaussian(chi);

  XMLFileWriter xml("t_overbu.xml");
  push(xml,"t_overbu");

  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
  write(xml,"psi", psi);
  write(xml,"chi", chi);

  Real OverMass = 1.5;
  Real m_q = 0.2;

  OverlapBULinOp over(u,OverMass,m_q);

  Real RsdCG = 1.0e-5;
  int  MaxCG = 1000;
  int  n_count;

  // Solve   D^dag.D*psi = D^dag*chi
  LatticeFermion chi_tmp;
  over(chi_tmp, chi, MINUS);
  InvCG2(over, chi_tmp, psi, RsdCG, MaxCG, n_count);

  // Check solution
  LatticeFermion  tmp;
  over(tmp, psi, PLUS);
  Double solnorm = sqrt(norm2(tmp-chi));
  Double ratio = solnorm / sqrt(norm2(chi));
  QDPIO::cout << "|solution| / |source| = " << ratio << endl;

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
