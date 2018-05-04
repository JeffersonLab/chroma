
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

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;

  int iter = 10000;
  QDPIO::cout << "Iters is " << iter << std::endl;

  //! Create a linear operator
  QDPIO::cout << "Constructing WilsonDslash" << std::endl;

  // WilsonDslash class can be optimised
  WilsonDslash D(u);

  QDPIO::cout << "Done" << std::endl;

  int i;

  int isign, cb, loop;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      clock_t myt1;
      clock_t myt2;
      double mydt;
      
      QDPIO::cout << "Applying D" << std::endl;
      
      myt1=clock();
      for(i=0; i < iter; i++) { 
	D.apply(chi, psi, (isign == 1 ? PLUS : MINUS), cb);
      }
      myt2=clock();
      
      mydt=(double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
      mydt=1.0e6*mydt/((double)(iter*(Layout::vol()/2)));
      
      QDPIO::cout << "cb = " << cb << " isign = " << isign << std::endl;
      QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		  << " (" <<  (double)(1392.0f/mydt) << ") Mflops " << std::endl;
    }
  }
  

  //! Create and try a more sophisticated operator
  /* Real Kappa = 0.1;
  PreconditionedWilson  M(u,Kappa);
  LatticeFermion eta;
  M(eta, psi, PLUS);

  write(nml,"eta", eta);
  */

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
