// $Id: t_follana_io_s.cc,v 1.4 2003-09-11 16:10:46 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "io/follana_io.h"

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

  XMLReader xml_in("input.xml");

  // XML File will look like this
  // 
  // <? xml version="1.0" ?>
  // <test>
  //    <nrow>X Y Z T</nrow>
  //    <infile>filename_of_propagator</infile>
  // </test>
  //
  multi1d<int> nrow(Nd);
  string filename;

  read(xml_in, "/test/nrow", nrow);
  read(xml_in, "/test/infile", filename);



  Layout::setLattSize(nrow);
  Layout::create();

  if ( Layout::primaryNode() ) { 
    cout << "Lattice: Lx = " << nrow[0] << " Ly = " << nrow[1] << " Lz = " << nrow[2]
	 << " Lt =" << nrow[3] << endl;

    cout << "Reading Propagator from file " << filename << endl;
  }

  // Try and read the propagator;

  LatticePropagator qprop= zero;
  readQpropFollana((char *)filename.c_str(), qprop);

  LatticeComplex corr_fn = trace(adj(qprop)*qprop);

  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(3));

  /* Project on zero momentum: Do a slice-wise sum. */
  multi1d<DComplex> hsum(nrow[3]);
  
  hsum = sumMulti(corr_fn, timeslice);

  XMLFileWriter xml_out("output.xml");

  push(xml_out, "follanaIO");
  Write(xml_out, hsum);
  pop(xml_out);

  // Time to bolt
  QDP_finalize();
  exit(0);

}
