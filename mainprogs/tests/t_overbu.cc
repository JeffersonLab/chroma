// $Id: t_overbu.cc,v 1.4 2003-09-11 00:46:04 edwards Exp $

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

  NmlWriter nml("t_overbu.nml");
  Write(nml,Nd);
  Write(nml,Nc);
  Write(nml,Ns);
  Write(nml,nrow);
  Write(nml,psi);
  Write(nml,chi);

  Real OverMass = 1.5;
  Real m_q = 0.2;

  OverlapBULinOp over(u,OverMass,m_q);

  Real RsdCG = 1.0e-5;
  int  MaxCG = 1000;
  int  n_count;

  // Solve   D^dag.D*psi = D^dag*chi
  LatticeFermion chi_tmp = over(chi, MINUS);
  InvCG2(over, chi_tmp, psi, RsdCG, MaxCG, n_count);

  // Check solution
  Double solnorm = sqrt(norm2(over(psi,PLUS)-chi));
  Double ratio = solnorm / sqrt(norm2(chi));
  cout << "|solution| / |source| = " << ratio << endl;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
