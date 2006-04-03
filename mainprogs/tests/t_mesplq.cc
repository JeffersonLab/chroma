// $Id: t_mesplq.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

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

  XMLFileWriter xml("t_mesplq.xml");
  push(xml, "t_mesplq");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;

#if 1
  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);
#else
  {
    XMLReader gauge_xml;
    readSzin(gauge_xml, u, string("CFGIN"));
    xml << gauge_xml;
  }

  unitarityCheck(u);
#endif

  // Try out the plaquette routine
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  // Test polyakov routine
  multi1d<DComplex> pollp(Nd);
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u, pollp[mu], mu);

  // Write out the results
  push(xml,"Observables");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  write(xml,"pollp", pollp);
  pop(xml);


  // Test gauge invariance
  rgauge(u);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  for(int mu = 0; mu < Nd; ++mu)
    polylp(u, pollp[mu], mu);

  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  push(xml,"Observables_after_gt");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  write(xml,"pollp", pollp);
  pop(xml);

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
