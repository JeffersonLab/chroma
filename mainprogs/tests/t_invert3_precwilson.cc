// $Id: t_invert3_precwilson.cc,v 1.1 2004-03-22 15:29:30 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "actions/ferm/invert/invcg2_timing_hacks.h"

using namespace QDP;
using namespace std;


enum GaugeStartType { COLD_START=0, HOT_START };
struct Params_t { 
  multi1d<int> nrow;
  multi1d<int> boundary;
  GaugeStartType gauge_start_type;
};

  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  Params_t params;

  // Read params
  XMLReader reader("DATA");
  string stype;
  try { 
    read(reader, "/t_invert/params/nrow", params.nrow);
    read(reader, "/t_invert/params/gauge_start_type", stype);
  }
  catch(const string &error) { 
    QDPIO::cerr << "Error : " << error << endl;
    throw;
  }
  reader.close();

  if( stype == "HOT_START" ) { 
    params.gauge_start_type = HOT_START;
  }
  else if( stype == "COLD_START" ) { 
    params.gauge_start_type = COLD_START;
  }

  QDPIO::cout << "Gauge start type " ;
  switch (params.gauge_start_type) { 
  case HOT_START:
    QDPIO::cout << "hot start" << endl;
    break;
  case COLD_START:
    QDPIO::cout << "cold start" << endl;
    break;
  default:
    QDPIO::cout << endl;
    QDPIO::cerr << "Unknown gauge start type " << endl;
  }

  params.boundary.resize(4);
  params.boundary[0] = 1;
  params.boundary[1] = 1;
  params.boundary[2] = 1;
  params.boundary[3] = 1;


  // Setup the lattice
  Layout::setLattSize(params.nrow);
  Layout::create();

  XMLFileWriter xml("XMLDAT");
  push(xml,"t_invert");
  push(xml,"params");
  write(xml, "nrow", params.nrow);
  write(xml, "boundary", params.boundary);
  write(xml, "gauge_start_type", stype);
  pop(xml); // Params

  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(params.boundary));
  
  // The Gauge Field
  multi1d<LatticeColorMatrix> u(Nd);
  
  switch ((GaugeStartType)params.gauge_start_type) { 
  case COLD_START:
    for(int j = 0; j < Nd; j++) { 
      u(j) = Real(1);
    }
    break;
  case HOT_START:
    // Hot (disordered) start
    for(int j=0; j < Nd; j++) { 
      random(u(j));
      reunit(u(j));
    }
    break;
  default:
    ostringstream startup_error;
    startup_error << "Unknown start type " << params.gauge_start_type <<endl;
    throw startup_error.str();
  }


  // Measure the plaquette on the gauge
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml, "plaquette");
  write(xml, "w_plaq", w_plaq);
  write(xml, "s_plaq", s_plaq);
  write(xml, "t_plaq", t_plaq);
  write(xml, "link", link);
  pop(xml);

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.
  EvenOddPrecWilsonFermAct  S_w(fbc, Real(0));

  Handle<const ConnectState> connect_state(S_w.createState(u));


  // Make me a linop (this callls the initialise function)
  Handle<const LinearOperator<LatticeFermion> > D_op(S_w.linOp(connect_state));

  
  LatticeFermion chi;
  LatticeFermion psi;
  StopWatch swatch;
  Double mydt;

#if 0
  LatticeFermion bchi,bpsi;
  gaussian(bchi);
  gaussian(bpsi);

  (*D_op)(bchi,bpsi,PLUS);
  
  exit(0);
#endif

  int iter;
  for(iter=1; ; iter <<= 1)
  {
    psi = zero;
    QDPIO::cout << "Let 0 action inverter iterate "<< iter << " times" << endl;

    gaussian(chi);
    swatch.reset();
    swatch.start();

    InvCG2_timing_hacks(*D_op, 
			chi, 
			psi,
			Real(1.0e-7),
			100000,
			iter);

    swatch.stop();
                                                                                
    mydt=Double(swatch.getTimeInSeconds());

    Internal::globalSum(mydt);

    mydt /= Double(Layout::numNodes());
                                                                                
    QDPIO::cout << "Time was " << mydt << " seconds" << endl;

    if ( toBool(mydt > Double(1) ) )
      break;
  }


  // Snarfed from SZIN
  int  N_dslash = 1320;
  int  N_mpsi   = 2*24 + 2*N_dslash;
  int  Nflops_c = (24 + 2*N_mpsi) + (48);     
  int Nflops_s = (2*N_mpsi + (2*48+2*24));   
  Double Nflops;

  multi1d<Double> mflops(10);
  multi1d<Double> mydt_a(10);
  
  for (int j=0; j < 10; ++j)
  {
    psi = zero;

    swatch.reset();
    swatch.start();

                                                                                    InvCG2_timing_hacks(*D_op, 
			chi, 
			psi,
			Real(1.0e-7),
			100000,
			iter);

    swatch.stop();
										    mydt=Double(swatch.getTimeInSeconds());

    Internal::globalSum(mydt);

    mydt /= Double(Layout::numNodes());
                                                                                
    mydt_a[j] = Double(1.0e6)*mydt/(Double(Layout::sitesOnNode())/Double(2));
                                                                                
    // Flop count for inverter 
    Nflops   = Double(Nflops_c) + Double(iter*Nflops_s);
    mflops[j] = Nflops / mydt_a[j];
  }

  mydt=1.0e6f*mydt/((Double)(Layout::sitesOnNode())/Double(2));

  push(xml, "Time");
  write(xml, "iter", iter);
  write(xml, "N_dslash", N_dslash);
  write(xml, "N_mpsi", N_mpsi);
  write(xml, "Nflops", Nflops);
  write(xml, "mydt_a", mydt_a);
  write(xml, "mflops_a", mflops);
  pop(xml);
  
  pop(xml);
  xml.close();

  QDP_finalize();
    
  exit(0);
}
