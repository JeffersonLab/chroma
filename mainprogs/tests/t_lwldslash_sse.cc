// $Id: t_lwldslash_sse.cc,v 1.7 2003-09-13 11:03:09 bjoo Exp $

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

  // Setup the layout
  XMLReader xml_in("input.xml");
  multi1d<int> nrow(Nd);
 
  read(xml_in, "/param/nrow", nrow);
  xml_in.close();

  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_lwldslash.xml");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi, chi2;
  gaussian(psi);
  chi = zero;

  int iter = 1000;
  cout << "Iters is " << iter << endl;

  //! Create a linear operator
  if( Layout::primaryNode() ) { 
    cout << "Constructing WilsonDslash" << endl;
  }

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
	cout << "The time per lattice point is "<< mydt 
	     << " micro sec (" <<  (double)(1392.0f/mydt) << ") Mflops " << endl;
	
	
      }
    }
  }
  
  //! Create a linear operator
  if( Layout::primaryNode() ) { 
    cout << "Constructing SSEWilsonDslash" << endl;
  }

  SSEWilsonDslash D_sse(u);

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
	chi = D_sse.apply(psi, (isign == 1 ? PLUS : MINUS ) , cb);
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
  
  // Correctness and consistency test
  // chi = D.apply(psi, PLUS, 0);
  gaussian(psi);
  Double n2;
  // Fill chi1 and make chi2 equal to it.
  LatticeFermion chi3;
  gaussian(chi3);
  
  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 0);
  chi2 = D.apply(psi, PLUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, +, 0) - D(psi, +, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 0);
  chi2 = D.apply(psi, PLUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, +, 1) - D(psi, +, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 1);
  chi2 = D.apply(psi, PLUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, +, 0) - D(psi, +, 1) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 1);
  chi2 = D.apply(psi, PLUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, +, 1) - D(psi, +, 1) || = " << n2 << endl;
  } 


  // MINUS

  chi = chi3; 
  chi2 = chi3;
  

  chi = D.apply(psi, MINUS, 0);
  chi2 = D.apply(psi, MINUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, -, 0) - D(psi, -, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 0);
  chi2 = D.apply(psi, MINUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, -, 1) - D(psi, -, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 1);
  chi2 = D.apply(psi, MINUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, -, 0) - D(psi, -, 1) || = " << n2 << endl;
  } 


  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 1);
  chi2 = D.apply(psi, MINUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D(psi, -, 1) - D(psi, -, 1) || = " << n2 << endl;
  } 

  // PLUS
  
  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, PLUS, 0);
  chi2 = D_sse.apply(psi, PLUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 0) - D_sse(psi, +, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, PLUS, 0);
  chi2 = D_sse.apply(psi, PLUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 1) - D_sse(psi, +, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, PLUS, 1);
  chi2 = D_sse.apply(psi, PLUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 0) - D_sse(psi, +, 1) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, PLUS, 1);
  chi2 = D_sse.apply(psi, PLUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 1) - D_sse(psi, +, 1) || = " << n2 << endl;
  } 


  // MINUS

  
  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, MINUS, 0);
  chi2 = D_sse.apply(psi, MINUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 0) - D_sse(psi, -, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, MINUS, 0);
  chi2 = D_sse.apply(psi, MINUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 1) - D_sse(psi, -, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, MINUS, 1);
  chi2 = D_sse.apply(psi, MINUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 0) - D_sse(psi, -, 1) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D_sse.apply(psi, MINUS, 1);
  chi2 = D_sse.apply(psi, MINUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 1) - D_sse(psi, -, 1) || = " << n2 << endl;
  } 


  // Mixed
  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 0);
  chi2 = D_sse.apply(psi, PLUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 0) - D(psi, +, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 0);
  chi2 = D_sse.apply(psi, PLUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 1) - D(psi, +, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 1);
  chi2 = D_sse.apply(psi, PLUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 0) - D(psi, +, 1) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, PLUS, 1);
  chi2 = D_sse.apply(psi, PLUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, +, 1) - D(psi, +, 1) || = " << n2 << endl;
  } 


  // MINUS

  
  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 0);
  chi2 = D_sse.apply(psi, MINUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 0) - D(psi, -, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 0);
  chi2 = D_sse.apply(psi, MINUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 1) - D(psi, -, 0) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 1);
  chi2 = D_sse.apply(psi, MINUS, 0); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 0) - D(psi, -, 1) || = " << n2 << endl;
  } 

  chi = chi3; 
  chi2 = chi3;

  chi = D.apply(psi, MINUS, 1);
  chi2 = D_sse.apply(psi, MINUS, 1); 
  n2 = norm2( chi2 - chi );

  if( Layout::primaryNode() ) { 
    cout << "|| D_sse(psi, -, 1) - D(psi, -, 1) || = " << n2 << endl;
  } 

  chi =zero;
  psi= zero;
  chi = D.apply(psi, PLUS, 0);

  NmlWriter nml("dump.nml");
  push(nml, "vectors");
  Write(nml, psi);
  Write(nml, chi);
  pop(nml);

  chi = zero;
  psi = zero;
  chi = D_sse.apply(psi, PLUS, 0);
  NmlWriter nml2("dump_sse.nml");
  push(nml2, "vectors");
  Write(nml2, psi);
  Write(nml2, chi);
  pop(nml2);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
