// $Id: t_lwldslash_array.cc,v 3.3 2007-03-07 00:03:01 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);


  int N5;

  // Lattice Size
  multi1d<int> nrow(Nd);
  try {
    // Read parameters
    XMLReader xml_in(Chroma::getXMLInputFileName());
    read(xml_in, "/param/nrow", nrow);
    read(xml_in, "/param/N5", N5);
    xml_in.close(); 
  }
  catch( const std::string&e ) { 
    QDPIO::cerr << "Caught Exception while reading XML " << e << endl;
    QDP_abort(1);
  }

  // Setup the layout
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml(Chroma::getXMLOutputFileName());
  push(xml,"t_lwldslash_array");
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
  QDPWilsonDslashOpt D(state);

  QDPIO::cout << "Done" << endl;


  push(xml,"Unoptimized_test");

  int isign, cb, loop, iter=1;

  bool first = true;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      double mydt; 
      QDP::StopWatch swatch;

      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying naive D " << iter << " times" << endl;

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
	
      QDPIO::cout << "Applying naive D for timings" << endl;
     
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

      push(xml,"Unopt_test");
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

      double mydt= 2.0;
      QDP::StopWatch swatch;
 
      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying D_opt " << iter << " times" << endl;

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

      QDPIO::cout << "Applying D_opt for timings" << endl;
     
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

      push(xml,"OPT_test");
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

      push(xml,"OPT_correctness_test");
      write(xml,"isign", isign);
      write(xml,"cb", cb);
      write(xml,"norm2_diff",n2);
      pop(xml);
    }
  }


  // Now do the vector case:
  multi1d<LatticeFermion> chis1(N5);
  multi1d<LatticeFermion> chis2(N5);
  multi1d<LatticeFermion> chis3(N5);
  multi1d<LatticeFermion> psis(N5);

  // Naive Dslash
  QDPIO::cout << "Consturcting Naive 5D Dslash, N5=" << N5 << endl;
  QDPWilsonDslashArrayOpt D5(state,N5);
  QDPIO::cout << "Done" << endl;


  push(xml,"Unoptimized_array_test");
  for(int k=0; k < N5; k++) { 
    gaussian(psis[k]);
    chis1[k]=zero;
    chis2[k]=zero;
    chis3[k]=zero;
  }

  iter=1;
  first = true;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      QDP::StopWatch swatch;
      double mydt=2.0;
     
      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying naive D5 " << iter << " times" << endl;

	  swatch.reset();
          swatch.start();
	  for(int i=iter; i-- > 0; ) {
	    D5.apply(chis1, psis, (isign == 1 ? PLUS : MINUS), cb);
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
	
      QDPIO::cout << "Applying naive D5 for timings" << endl;
     
      swatch.reset();
      swatch.start(); 
      for(int i=iter; i-- > 0; ) {
	D5.apply(chis1, psis, (isign == 1 ? PLUS : MINUS), cb);
      }
      swatch.stop();
      
      mydt=swatch.getTimeInSeconds();
      mydt=1.0e6*mydt/double(iter*(Layout::sitesOnNode()/2));
      QDPInternal::globalSum(mydt);
      mydt /= Layout::numNodes();
 
      float mflops = float(1320.0f*N5/mydt);
      QDPIO::cout << "cb = " << cb << " isign = " << isign << endl;
      QDPIO::cout << "The time per lattice point is "<< mydt 
		  << " micro sec (" <<  mflops << ") Mflops " << endl;

      push(xml,"Unopt_array_test");
      write(xml,"cb",cb);
      write(xml,"isign",isign);
      write(xml,"mflops",mflops);
      pop(xml);
    }
  }
  
  pop(xml);

  //! Naive loop
  QDPIO::cout << "Constructing (possibly optimized) WilsonDslash to do vector operation with a loop" << endl;
  WilsonDslash D_opt_loop(state);
  QDPIO::cout << "Done" << endl;

  push(xml,"Optimized_loop_test");

  first = true;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      QDP::StopWatch swatch;
      double mydt= 2.0;
      
      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying D_opt_loop " << iter << " times" << endl;

	  swatch.reset();
	  swatch.start();
	  for(int i=iter; i-- > 0; ) {
	    for(int loop=0; loop < N5; loop++) { 
	      D_opt_loop.apply(chis1[loop], psis[loop], (isign == 1 ? PLUS : MINUS ) , cb); // NOTE: for timings throw away return value
	    }
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

      QDPIO::cout << "Applying D_opt_loop for timings" << endl;
      
      swatch.reset();
      swatch.start();
      for(int i=iter; i-- > 0; ) {
	for(int loop=0; loop<N5; loop++) { 
	  D_opt_loop.apply(chis1[loop], psis[loop], (isign == 1 ? PLUS : MINUS ) , cb); // NOTE: for timings throw away return value
	}
      }
      swatch.stop();
      
      mydt=swatch.getTimeInSeconds();
      mydt=1.0e6*mydt/double(iter*(Layout::sitesOnNode()/2));
      QDPInternal::globalSum(mydt);
      mydt /= Layout::numNodes();
 
      float mflops = float(1320.0f*N5/mydt);
      QDPIO::cout << "cb = " << cb << " isign = " << isign << endl;
      QDPIO::cout << "After " << iter << " calls, the time per lattice point is "<< mydt 
		  << " micro sec (" <<  mflops << ") Mflops " << endl;

      push(xml,"OPT_loop_test");
      write(xml,"cb",cb);
      write(xml,"isign",isign);
      write(xml,"mflops",mflops);
      pop(xml);
    }
  }

  pop(xml);

  //! Create a linear operator
  QDPIO::cout << "Constructing (possibly optimized) WilsonDslashArray N5="<< N5 << endl;

  WilsonDslashArray D5_opt(state, N5);

  QDPIO::cout << "Done" << endl;

  push(xml,"Optimized_array_test");

  first = true;
  for(isign = 1; isign >= -1; isign -= 2) {
    for(cb = 0; cb < 2; ++cb) { 

      double mydt= 2.0;
      QDP::StopWatch swatch;
 
      if (first) 
      {
	for(iter=1; ; iter <<= 1)
	{
	  QDPIO::cout << "Applying D5_opt " << iter << " times" << endl;

	  swatch.reset();
	  swatch.start();
	  for(int i=iter; i-- > 0; ) {
	    D5_opt.apply(chis1, psis, (isign == 1 ? PLUS : MINUS ) , cb); // NOTE: for timings throw away return value
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

      QDPIO::cout << "Applying D5_opt for timings" << endl;
      
      swatch.reset();
      swatch.start(); 
      for(int i=iter; i-- > 0; ) {
	D5_opt.apply(chis1, psis, (isign == 1 ? PLUS : MINUS ) , cb); // NOTE: for timings throw away return value
      }
      swatch.stop();
      
      mydt=swatch.getTimeInSeconds();
      mydt=1.0e6*mydt/double(iter*(Layout::sitesOnNode()/2));
      QDPInternal::globalSum(mydt);
      mydt /= Layout::numNodes();
 
      float mflops = float(1320.0f*N5/mydt);
      QDPIO::cout << "cb = " << cb << " isign = " << isign << endl;
      QDPIO::cout << "After " << iter << " calls, the time per lattice point is "<< mydt 
		  << " micro sec (" <<  mflops << ") Mflops " << endl;

      push(xml,"Opt_array_test");
      write(xml,"cb",cb);
      write(xml,"isign",isign);
      write(xml,"mflops",mflops);
      pop(xml);
    }
  }

  pop(xml);

  for(cb = 0; cb < 2; cb++) { 
    for(isign = 1; isign >= -1; isign -= 2) { 
      for(int k=0; k < N5; k++) { 
	chis1[k] = zero;
	chis2[k] = zero;
      }

      D5.apply(chis1, psis, (isign > 0 ? PLUS : MINUS), cb);
      D5.apply(chis2, psis, (isign > 0 ? PLUS : MINUS), cb);
      
      // 5D Norm
      n2 = Double(0);
      for(int i=0; i < N5; i++) { 
	n2 += norm2( chis2[i] - chis1[i] );
      }

      QDPIO::cout << "Paranoia test: || D5(psi, "
		  << (isign > 0 ? "+, " : "-, ") <<  cb 
		  << ") - D5(psi, " 
		  << (isign > 0 ? "+, " : "-, ") <<  cb << " ) ||  = " << n2 
		  << endl;
    }
  }

  for(cb = 0; cb < 2; cb++) { 
    for(isign = 1; isign >= -1; isign -= 2) { 

      for(int k=0; k < N5; k++) { 
	chis1[k] = zero;
	chis2[k] = zero;
      }

      D5.apply(chis1, psis, (isign > 0 ? PLUS : MINUS), cb);
      D5_opt.apply(chis2, psis, (isign > 0 ? PLUS : MINUS), cb);

      n2=Double(0);
      for(int i=0; i < N5; i++) {
	n2 += norm2( chis2[i] - chis1[i] );
      }
      QDPIO::cout << "OPT test: || D5(psi, "
		  << (isign > 0 ? "+, " : "-, ") <<  cb 
		  << ") - D5_opt(psi, " 
		  << (isign > 0 ? "+, " : "-, ") <<  cb << " ) ||  = " << n2 
		  << endl;

      push(xml,"Opt_array_correctness_test");
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
