// $Id: t_follana_io_s.cc,v 1.2 2003-09-11 12:31:20 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "follana_io.h"

using namespace QDP;



//! Function object used for constructing the time-slice set
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}
                                                                                
  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}
                                                                                
  int dir_decay;
                                                                                
private:
  TimeSliceFunc() {}  // hide default constructor
};



int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {16,16,16,32};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements


  Layout::setLattSize(nrow);
  Layout::create();

  // Try and read the propagator;

  LatticePropagator qprop;
  readQpropFollana("./prop.0", qprop);

  LatticeComplex corr_fn = trace(adj(qprop)*qprop);

  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(3));

  /* Project on zero momentum: Do a slice-wise sum. */
  multi1d<DComplex> hsum(nrow[3]);
  
  
  hsum = sumMulti(corr_fn, timeslice);


  NmlWriter nml("prop.nml");

  // push(nml, "propagator");
  //Write(nml, qprop);
  // pop(nml);
  push(nml, "goldstone");
  Write(nml,hsum);
  pop(nml);

  // Time to bolt
  QDP_finalize();
}
