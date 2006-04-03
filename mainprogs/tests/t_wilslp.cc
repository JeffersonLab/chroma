// $Id: t_wilslp.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $
//
// Test driver for wilslp
// default test uses the usual NERSC 4^4 config.
// 


#include <iostream>
#include <cstdio>
#include <time.h>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
    const int foo[] = {4,4,4,4}; 
//  const int foo[] = {16,16,16,32}; 
//  const int foo[] = {20,20,20,64}; 

  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_wilslp.xml");
  push(xml, "t_wilslp");

  push(xml,"lattice");
  write(xml,"Nd",Nd);
  write(xml,"Nc",Nc);
  write(xml,"nrow",nrow);
  pop(xml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);
  Double w_plaq, s_plaq, t_plaq, link;

/*
  QDPIO::cerr << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);
*/

 XMLReader gauge_xml;
 string gaugeFilename;

 gaugeFilename="../../tests/t_asqtad_prop/t_nersc.cfg";
// gaugeFilename="/raidz/aci/configs/getlat.1200";
// gaugeFilename="/raidz/aci/configs/u_MILC_2064f21b681m030m050.246";

  QDPIO::cout << "Reading config: " << gaugeFilename << endl;
   readArchiv(gauge_xml, u, gaugeFilename);

  // Reunitarize the gauge field
  for(int m=0; m < u.size(); ++m)
    reunit(u[m]);

// Check gauge invariance here if required
//    rgauge(u);
//    QDPIO::cout << "Random gauge transformation for testing " << endl;

  // Try out the plaquette routine
  QDPIO::cout << "Measure the plaquette and link " << endl;

  StopWatch swatch;
  swatch.start();

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  swatch.stop();
  double time_in_sec  = swatch.getTimeInSeconds();

  QDPIO::cout << "w_plaq = " << w_plaq << endl;
  QDPIO::cout << "s_plaq = " << s_plaq << endl;
  QDPIO::cout << "t_plaq = " << t_plaq << endl;
  QDPIO::cout << "link = " << link << endl;

  QDPIO::cout << "Measurement took " << time_in_sec 
	<< " secs" << endl;

  // Write out the results
  push(xml,"observables");
  write(xml,"w_plaq",w_plaq);
  write(xml,"w_link",link);


  // Now try calling the Wilson loop routine
  QDPIO::cout << "Measure Wilson loops" << endl;
  int j_decay = Nd - 1;
  int t_dir = j_decay;

  int kind    = 6;
// kind -- binary-combined YES/NO [1/0] of the three loop options     
//  e.g. kind = 2 gives planar t-like
//              6       planar + off-axis: sqrt(2), sqrt(5), sqrt(3) 


//  Real sm_fact = 2.5;   // typical parameter
//  int n_smear = 0;     // number of smearing hits

//  int BlkMax = 100;    // max iterations in max-ing trace
// Real BlkAccu = 1.0e-5;  // accuracy of max-ing
  
//  QDPIO::cout << "sm_fact: " << sm_fact << endl;
//  QDPIO::cout << "n_smear: " << n_smear << endl;

  swatch.reset();
  swatch.start();

  wilslp(u, j_decay, t_dir, kind,
	 xml, "Fuzzed_Wilson_Loops");

  swatch.stop();
  time_in_sec  = swatch.getTimeInSeconds();
  QDPIO::cout << "Loop measurements complete" << endl;
  QDPIO::cout << "Measurement took " << time_in_sec 
	<< " secs" << endl;

  pop(xml);
  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
