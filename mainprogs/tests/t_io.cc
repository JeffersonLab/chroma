// $Id: t_io.cc,v 1.8 2004-02-11 12:51:35 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  LatticeReal a;
  Double d = 17;
  random(a);

  NmlWriter tonml("cat");
  write(tonml,"a", a);
  tonml.close();

  BinaryWriter tobinary("dog");
  write(tobinary, a);
  write(tobinary, d);
  tobinary.close();

  float x;
  QDPIO::cout << "Read some data from file input\n";
  TextReader from("input");
  from >> x;
  from.close();

  QDPIO::cout << "you entered :" << x << ":\n";
  
  // Time to bolt
  QDP_finalize();

  exit(0);
}
