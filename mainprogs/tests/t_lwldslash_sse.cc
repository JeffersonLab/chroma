// $Id: t_lwldslash_sse.cc,v 3.1 2006-06-13 20:15:25 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Read parameters
  XMLReader xml_in("input.xml");

  // Lattice Size
  multi1d<int> nrow(Nd);
  read(xml_in, "/param/nrow", nrow);

  xml_in.close();

  // Setup the layout
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_lwldslash.xml");

  push(xml,"t_lwldslash");

  proginfo(xml);    // Print out basic program info

  // Make up a random gauge field.
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);


  // Make up a gaussian source and a zero result vector
  LatticeFermion psi, chi, chi2;
  gaussian(psi);
  chi = zero;

  //! Create a linear operator
  QDPIO::cout << "Constructing naive QDPWilsonDslash" << endl;

  Handle< FermState<LatticeFermion,
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> > > state(new PeriodicFermState<LatticeFermion,
					  multi1d<LatticeColorMatrix>,
					  multi1d<LatticeColorMatrix> >(u));

  // Naive Dslash
  QDPWilsonDslash D(state);

  QDPIO::cout << "Done" << endl;


  push(xml,"Unoptimized_test");

  int isign, cb, loop, iter=1;
  bool first = true;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      QDP::StopWatch swatch;
      double mydt;
 
      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying D " << iter << " times" << endl;

	  swatch.reset();
  	  swatch.start();
	  for(int i=iter; i-- > 0; ) {
	    D.apply(chi, psi, (isign == 1 ? PLUS : MINUS), cb);
	  }
	  swatch.stop();
	  mydt=swatch.getTimeInSeconds();

	  QDPInternal::globalSum(mydt);
	  mydt /= Layout::numNodes();

	  if (mydt > 1) {
	    first = false;
	    break;
	  }
	}
      }
	
      QDPIO::cout << "Applying D for timings" << endl;
    
      swatch.reset();
      swatch.start(); 
      for(int i=iter; i-- > 0; ) {
	D.apply(chi, psi, (isign == 1 ? PLUS : MINUS), cb);
      }
      swatch.stop();
      
      mydt=swatch.getTimeInSeconds();
      mydt=1.0e6*mydt/double(iter*(Layout::sitesOnNode()/2));
      QDPInternal::globalSum(mydt);
      mydt /= Layout::numNodes();
 
      float mflops = float(1320.0f/mydt);
      QDPIO::cout << "cb = " << cb << " isign = " << isign << endl;
      QDPIO::cout << "The time per lattice point is "<< mydt 
		  << " micro sec (" <<  mflops << ") Mflops " << endl;

      push(xml,"test");
      write(xml,"cb",cb);
      write(xml,"isign",isign);
      write(xml,"mflops",mflops);
      pop(xml);
    }
  }
  
  pop(xml);

  //! Create a linear operator
  QDPIO::cout << "Constructing (possibly optimized) WilsonDslash" << endl;

  WilsonDslash D_opt(state);

  QDPIO::cout << "Done" << endl;

  push(xml,"Optimized_test");

  first = true;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      double mydt;
      QDP::StopWatch swatch;
 
      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying D " << iter << " times" << endl;
	
	  swatch.reset();
	  swatch.start();

	  for(int i=iter; i-- > 0; ) {
	    D_opt.apply(chi, psi, (isign == 1 ? PLUS : MINUS ) , cb); // NOTE: for timings throw away return value
	  }
	  swatch.stop();


	  mydt=swatch.getTimeInSeconds();
	  QDPInternal::globalSum(mydt);
	  mydt /= Layout::numNodes();

	  if (mydt > 1) {
	    first = false;
	    break;
	  }
	}
      }

      QDPIO::cout << "Applying D for timings" << endl;
     
      swatch.reset();
      swatch.start(); 
      for(int i=iter; i-- > 0; ) {
	D_opt.apply(chi, psi, (isign == 1 ? PLUS : MINUS ) , cb); // NOTE: for timings throw away return value
      }
      swatch.stop(); 
      mydt=swatch.getTimeInSeconds();
      mydt=1.0e6*mydt/double(iter*(Layout::sitesOnNode()/2));
      QDPInternal::globalSum(mydt);
      mydt /= Layout::numNodes();
 
      float mflops = float(1320.0f/mydt);
      QDPIO::cout << "cb = " << cb << " isign = " << isign << endl;
      QDPIO::cout << "After " << iter << " calls, the time per lattice point is "<< mydt 
		  << " micro sec (" <<  mflops << ") Mflops " << endl;

      push(xml,"test");
      write(xml,"cb",cb);
      write(xml,"isign",isign);
      write(xml,"mflops",mflops);
      pop(xml);
    }
  }

  pop(xml);

  LatticeFermion chi3;
  Double n2;

  gaussian(chi3);
  gaussian(psi);
  for(cb = 0; cb < 2; cb++) { 
    for(isign = 1; isign >= -1; isign -= 2) { 

      chi = chi3;
      chi2 = chi3;
      D.apply(chi, psi, (isign > 0 ? PLUS : MINUS), cb);
      D.apply(chi2, psi, (isign > 0 ? PLUS : MINUS), cb);
      
      n2 = norm2( chi2 - chi );

      QDPIO::cout << "Paranoia test: || D(psi, "
		  << (isign > 0 ? "+, " : "-, ") <<  cb 
		  << ") - D(psi, " 
		  << (isign > 0 ? "+, " : "-, ") <<  cb << " ) ||  = " << n2 
		  << endl;
    }
  }

  gaussian(chi3);
  gaussian(psi);
  for(cb = 0; cb < 2; cb++) { 
    for(isign = 1; isign >= -1; isign -= 2) { 

      chi = chi3;
      chi2 = chi3;
      D.apply(chi, psi, (isign > 0 ? PLUS : MINUS), cb);
      D_opt.apply(chi2, psi, (isign > 0 ? PLUS : MINUS), cb);
      
      n2 = norm2( chi2 - chi );

      QDPIO::cout << "OPT test: || D(psi, "
		  << (isign > 0 ? "+, " : "-, ") <<  cb 
		  << ") - D_opt(psi, " 
		  << (isign > 0 ? "+, " : "-, ") <<  cb << " ) ||  = " << n2 
		  << endl;

      push(xml,"OPT test");
      write(xml,"isign", isign);
      write(xml,"cb", cb);
      write(xml,"norm2_diff",n2);
      pop(xml);
    }
  }
	  
  pop(xml);
  
  // Time to bolt
  Chroma::finalize();

  exit(0);
}
