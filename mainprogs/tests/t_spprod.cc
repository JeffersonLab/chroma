// $Id: t_spprod.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_spprod.xml");
  push(xml,"t_spprod");

  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);

  LatticeFermion fpsi;
  LatticeFermion fchi;
  LatticeFermion fa1;
  LatticePropagator ppsi;
  LatticePropagator pchi;
  LatticePropagator pa1;

  /* Test 1 */
  printf("Fill fermions with gaussians and spprod\n");
  gaussian(fpsi);
  gaussian(fchi);

  push(xml,"here_is_psi");
  write(xml, "fpsi", fpsi);
  pop(xml);

  /* fa1 = (gamma(n))*psi */
  for(int n = 0; n < Ns*Ns; ++n)
  {
    fa1 = zero;
    fa1[rb[0]] = Gamma(n) * fpsi;

    printf("print the fa1 fields in direction n= %d\n", n);
    push(xml,"Basis");
    write(xml, "n", n);
    pop(xml);
    push(xml,"fa1_is_spin");
    write(xml, "fa1", fa1);
    pop(xml);
  }

      
  /* Test 2 */
  printf("Fill fermions with gaussians and spprod/trace\n");
  gaussian(fpsi);
  gaussian(fchi);

  /* fa1 = (gamma(n))*psi */
  for(int n = 0; n < Ns*Ns; ++n)
  {
    DComplex dcsum = innerProduct(fchi, Gamma(n) * fpsi);

    printf("print the fa1 fields in direction n= %d\n", n);
    push(xml,"Fermion_inner_product");
    write(xml, "n", n);
    write(xml, "dcsum", dcsum);
    pop(xml);
  }

          
  /* Test 3 */
  printf("Fill propagators with gaussians and spprod/trace\n");
  gaussian(ppsi);
  gaussian(pchi);

  /* pa1 = (gamma(n))*psi */
  for(int n = 0; n < Ns*Ns; ++n)
  {
    DComplex dcsum = innerProduct(pchi, Gamma(n) * ppsi);

    printf("print the pa1 fields in direction n= %d\n", n);
    push(xml,"Propagator_inner_product");
    write(xml, "n", n);
    write(xml, "dcsum", dcsum);
    pop(xml);
  }

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
