// $Id: t_lower_tests.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $
//
//  This is a collection of simple 
//  tests to make sure that none of 
//  the underlying routines have been 
//  broken.
//
//
//

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter& xml = Chroma::getXMLOutputInstance();
  push(xml, "collection_of_unit_tests");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

  // create a hot gauge gauge configuration
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;

  QDPIO::cerr << "Creating gauge configuration\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

  // Try out the plaquette routine
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  // Write out the results
  push(xml,"gauge_observables");
  write(xml,"config_type", "hot");
  write(xml,"w_plaq", w_plaq);
  write(xml,"link", link);
  pop(xml);

  QDPIO::cerr << "Starting tests\n";

  // 
  //  regression test the displacement operator
  //  

  push(xml,"simple_test");
  write(xml,"testname", "displacement");


  int length = 3 ; 
  int dir = 0 ; 

  LatticePropagator  chiA ; 
  LatticePropagator  chiB ; 

  gaussian(chiA) ; 
  chiB = chiA ; 
  displacement(u,chiA,length,dir); 

  Double tmp_ = sum(real(trace(chiA*chiA)));
  write(xml,"sumorg", tmp_);

  Double tmp = sum(real(trace(chiA*chiB)));
  write(xml,"sum", tmp);

  pop(xml);   // end of single test

  // final pop
  pop(xml);

  xml.close(); 

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
