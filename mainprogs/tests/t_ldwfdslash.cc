// $Id: t_ldwfdslash.cc,v 1.1 2003-10-20 20:34:04 edwards Exp $

#include <iostream>
#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>

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

  XMLFileWriter xml("t_lwldslash.xml");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeDWFermion psi, chi;
  random(psi);
  chi = zero;

  int iter = 10;
  QDPIO::cout << "Iters is " << iter << endl;

  //! Create a linear operator
  QDPIO::cout << "Constructing DWDslash" << endl;

  // DWDslash class can be optimised
  Real WilsonMass = 1.5;
  DWDslash D(u,WilsonMass);

  QDPIO::cout << "Done" << endl;

  for(int isign = 1; isign >= -1; isign -= 2) 
  {
    for(int cb = 0; cb < 2; ++cb) 
    { 
      QDPIO::cout << "Applying D" << endl;
      
      clock_t myt1=clock();
      for(int i=0; i < iter; i++)
	D.apply(psi, (isign == 1 ? PLUS : MINUS), cb);  // throw away the result
      clock_t myt2=clock();
      
      double mydt = (double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
      mydt *= 1.0e6/((double)(iter*(Layout::vol()/2)));
      
      QDPIO::cout << "cb = " << cb << " isign = " << isign << endl;
      QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		  << " (" <<  (double)(Ls*1392.0f/mydt) << ") Mflops " << endl;
    }
  }
  
  // Time to bolt
  QDP_finalize();

  exit(0);
}
