// $Id: t_lwldslash_sse.cc,v 1.10 2003-09-17 15:15:00 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "primitives.h" // GTF: for PLUS

#include <lib.h>

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Read parameters
  XMLReader xml_in("input.xml");

  // Lattice Size
  multi1d<int> nrow(Nd);
  read(xml_in, "/param/nrow", nrow);

  // No of iters
  int iter;
  read(xml_in, "/param/iter", iter);

  xml_in.close();

  // Setup the layout
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_lwldslash.xml");

  // Make up a random gauge field.
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  // Make up a gaussian source and a zero result vector
  LatticeFermion psi, chi, chi2;
  gaussian(psi);
  chi = zero;

  if( Layout::primaryNode() ) { 
    cout << "Going to do " << iter << " Dslash Applications." << endl;
  }

  //! Create a linear operator
  if( Layout::primaryNode() ) { 
    cout << "Constructing naive QDPWilsonDslash" << endl;
  }

  // Naive Dslash
  QDPWilsonDslash D(u);

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
	cout << "The time per lattice point is "<< mydt 
	     << " micro sec (" <<  (double)(1392.0f/mydt) << ") Mflops " << endl;
	
	
      }
    }
  }
  
  //! Create a linear operator
  if( Layout::primaryNode() ) { 
    cout << "Constructing (possibly optimized) WilsonDslash" << endl;
  }

  WilsonDslash D_opt(u);

  if( Layout::primaryNode() ) { 
    cout << "Done" << endl;
  }

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
	chi = D_opt.apply(psi, (isign == 1 ? PLUS : MINUS ) , cb);
      }
      myt2=clock();
      
      mydt=(double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
      mydt=1.0e6*mydt/((double)(iter*(Layout::vol()/2)));
      
      if( Layout::primaryNode() ) { 
	cout << "cb = " << cb << " isign = " << isign << endl;
	cout << "The time per lattice point is "<< mydt 
	     << " micro sec (" <<  (double)(1392.0f/mydt) << ") Mflops " << endl;
	
	
      }
    }
  }


  LatticeFermion chi3;
  Double n2;

  gaussian(chi3);
  gaussian(psi);
  for(cb = 0; cb < 2; cb++) { 
    for(isign = 1; isign >= -1; isign -= 2) { 

      chi = chi3;
      chi2 = chi3;
      chi[ rb[cb] ] = D.apply(psi, (isign > 0 ? PLUS : MINUS), cb);
      chi2[ rb[cb] ] = D.apply(psi, (isign > 0 ? PLUS : MINUS), cb);
      
      n2 = norm2( chi2 - chi );

      if( Layout::primaryNode() ) { 
	cout << "Paranoia test: || D(psi, "
	     << (isign > 0 ? "+, " : "-, ") <<  cb 
	     << ") - D(psi, " 
	     << (isign > 0 ? "+, " : "-, ") <<  cb << " ) ||  = " << n2 
	     << endl;
      }
    }
  }

  gaussian(chi3);
  gaussian(psi);
  for(cb = 0; cb < 2; cb++) { 
    for(isign = 1; isign >= -1; isign -= 2) { 

      chi = chi3;
      chi2 = chi3;
      chi[ rb[ cb ] ] = D_opt.apply(psi, (isign > 0 ? PLUS : MINUS), cb);
      chi2[ rb[ cb ] ] = D.apply(psi, (isign > 0 ? PLUS : MINUS), cb);
      
      n2 = norm2( chi2 - chi );

      if( Layout::primaryNode() ) { 
	cout << "OPT test: || D(psi, "
	     << (isign > 0 ? "+, " : "-, ") <<  cb 
	     << ") - D_opt(psi, " 
	     << (isign > 0 ? "+, " : "-, ") <<  cb << " ) ||  = " << n2 
	     << endl;
      }
    }
  }
	  
      
  
  // Time to bolt
  QDP_finalize();

  exit(0);
}
