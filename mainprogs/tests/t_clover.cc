// t_su3.cc, 2004/11/16 velytsky
// velytski@csit.fsu.edu
// gluodynamics (su(3)) heatbath updating with measurements
#include <iostream>
#include <cstdio>
#include <cmath>
#include "chroma.h"

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeColorMatrix> Q;
typedef multi1d<LatticeColorMatrix> P;


int main(int argc, char *argv[]) 
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[]={4, 4, 4, 4};
  multi1d<int> nrow(Nd);
  nrow=foo;
  Layout::setLattSize(nrow);
  Layout::create();

  // Start up a weak field
  struct Cfg_t config = { CFG_TYPE_WEAK_FIELD, "dummy" };
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, config);
  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Need a ferm state 
  Handle< FermState<T,P,Q> > fs(new PeriodicFermState<T,P,Q>(u));

  // Need some clover params 
  CloverFermActParams params;
  params.clovCoeffR = Real(1.92);
  params.clovCoeffT = Real(0.57);
  params.Mass=Real(0);

  // Need a clover operator - reference
  QDPCloverTerm qdp_clov;
  qdp_clov.create(fs, params);



  LatticeFermion src, dest1, dest2;
  gaussian(src);
  dest1=zero;
  dest2=zero;

  int n_sec=30;
  unsigned long iters=1;
  double time;

  QDPIO::cout << "Calibrating QDPCloverTerm" << endl;
  do { 
    iters *= 2;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    for(unsigned long i=0; i < iters; i++) { 
      qdp_clov.apply(dest1, src, PLUS, 0);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();
    QDPIO::cout << "." << flush;
  }
  while(time < (double)n_sec );


  QDPIO::cout << endl;

  {
    QDPIO::cout << "Timing with " << iters << " iterations" << endl;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    for(unsigned long i=0; i < iters; i++) { 
      qdp_clov.apply(dest1, src, PLUS, 0);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    // Flopcount
    QDPIO::cout << "Layout volume=" << Layout::vol() << endl;
    QDPIO::cout << "FLOPS=" << qdp_clov.nFlops() << endl;
    QDPIO::cout << "iters=" << iters <<endl;
    QDPIO::cout << "time=" << time << endl;
    double flops=(double)qdp_clov.nFlops();
    flops *= (double)(Layout::sitesOnNode())/(double)2;
    flops *= (double)(iters);
    FlopCounter fc;
    fc.addFlops((unsigned long long)floor(flops));
    fc.report("QDPCloverTerm", time);

  }

  CloverTerm opt_clov;
  opt_clov.create(fs, params);

  iters = 1;
  QDPIO::cout << "Calibrating CloverTerm" << endl;
  do { 
    iters *= 2;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    for(unsigned long i=0; i < iters; i++) { 
      opt_clov.apply(dest1, src, PLUS, 0);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();
    QDPIO::cout << "." << flush;
  }
  while(time < (double)n_sec );


  QDPIO::cout << endl;

  {
    QDPIO::cout << "Timing with " << iters << " iterations" << endl;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    for(unsigned long i=0; i < iters; i++) { 
      opt_clov.apply(dest1, src, PLUS, 0);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    // Flopcount
    QDPIO::cout << "Layout volume=" << Layout::vol() << endl;
    QDPIO::cout << "FLOPS=" << qdp_clov.nFlops() << endl;
    QDPIO::cout << "iters=" << iters <<endl;
    QDPIO::cout << "time=" << time << endl;
    double flops=(double)qdp_clov.nFlops();
    flops *= (double)(Layout::sitesOnNode())/(double)2;
    flops *= (double)(iters);
    FlopCounter fc;
    fc.addFlops((unsigned long long)floor(flops));
    fc.report("CloverTerm", time);

  }


  

  iters=1;
  EvenOddPrecCloverLinOp D;
  D.create(fs, params);

  QDPIO::cout << "Calibrating Possibly optimized EOPREC Clov Op" << endl;
  do { 
    iters *= 2;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iters; i++) { 
      D(dest1, src, PLUS);
      D(dest1, src, MINUS);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();
    QDPIO::cout << "." << flush;
  }
  while(time < (double)n_sec );
  QDPIO::cout << endl;

  {
    QDPIO::cout << "Timing with " << iters << " iterations" << endl;
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iters; i++) { 
      D(dest1, src, PLUS);
      D(dest1, src, MINUS);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    FlopCounter fc;
    fc.reset();
    unsigned long long flops = (unsigned long long)2;
    flops *= (unsigned long long)iters;
    flops *= (unsigned long long)D.nFlops();

    fc.addFlops(flops);
    fc.report("EOPrecCloverOp", time);
  }


  
  Chroma::finalize();
  return 1;
}
