// $Id: t_lwldslash.cc,v 1.9 2003-09-16 13:38:37 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "primitives.h" // GTF: for PLUS


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

  LatticeFermion psi, chi;
  random(psi);
  chi = zero;

  int iter = 10000;
  cout << "Iters is " << iter << endl;

  //! Create a linear operator
  if( Layout::primaryNode() ) { 
    cout << "Constructing WilsonDslash" << endl;
  }

  // WilsonDslash class can be optimised
  WilsonDslash D(u);

  if( Layout::primaryNode() ) { 
    cout << "Done" << endl;
  }

  int i;

  int isign, cb, loop;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      clock_t myt1;
      clock_t myt2;
      double mydt;
      
      if( Layout::primaryNode() ) { 
	cout << "Applying D" << endl;
      }
      
      myt1=clock();
      for(i=0; i < iter; i++) { 
	chi = D.apply(psi, (isign == 1 ? PLUS : MINUS ) , cb);
      }
      myt2=clock();
      
      mydt=(double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
      mydt=1.0e6*mydt/((double)(iter*(Layout::vol()/2)));
      
      if( Layout::primaryNode() ) { 
	cout << "cb = " << cb << " isign = " << isign << endl;
	cout << "The time per lattice point is "<< mydt << " micro sec" 
	     << " (" <<  (double)(1392.0f/mydt) << ") Mflops " << endl;
	
	
      }
    }
  }
  

  //! Create and try a more sophisticated operator
  /* Real Kappa = 0.1;
  PreconditionedWilson  M(u,Kappa);
  LatticeFermion eta;
  eta = M(psi, PLUS);

  Write(nml,eta);
  */

  // Time to bolt
  QDP_finalize();

  exit(0);
}
