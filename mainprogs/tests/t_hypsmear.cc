// $Id: t_hypsmear.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_hypsmear.xml");
  push(xml,"t_hypsmear");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;

//  QDPIO::cout << "Reading test_purgaug.cfg1\n";
//  Seed seed_old;
//  readSzin(u, 0, "../test_purgaug.cfg1", seed_old);
  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

  // Try out the plaquette routine
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  // Write out the results
  push(xml,"observables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  pop(xml);

  // Now hyp smear
  multi1d<LatticeColorMatrix> u_hyp(Nd);
  Real BlkAccu = 1.0e-5;
  int BlkMax = 100;
  Real alpha1 = 0.7;
  Real alpha2 = 0.6;
  Real alpha3 = 0.3;
  Hyp_Smear(u, u_hyp, alpha1, alpha2, alpha3, BlkAccu, BlkMax);

  // Measure again
  MesPlq(u_hyp, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  // Write out the results
  push(xml,"HYP_observables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  pop(xml);

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
