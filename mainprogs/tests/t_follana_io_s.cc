// $Id: t_follana_io_s.cc,v 3.1 2007-02-22 21:11:50 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "io/follana_io_s.h"
#include "meas/hadron/mesphas_s.h"

using namespace Chroma;



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
  Chroma::initialize(&argc, &argv);

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

  QDPIO::cout << "Lattice: Lx = " << nrow[0] << " Ly = " << nrow[1] << " Lz = " << nrow[2]
	      << " Lt =" << nrow[3] << endl;

  QDPIO::cout << "Reading Propagator from file " << filename << endl;

  // Try and read the propagator;

  LatticePropagator qprop= zero;
  readQpropFollana((char *)filename.c_str(), qprop);

  LatticeComplex corr_fn = trace(adj(qprop)*qprop);

  Set timeslice;
  timeslice.make(TimeSliceFunc(3));

  /* Project on zero momentum: Do a slice-wise sum. */
  multi1d<DComplex> hsum(nrow[3]);
  
  hsum = sumMulti(corr_fn, timeslice);

  XMLFileWriter xml_out("output.xml");

  push(xml_out, "follanaIO");
  write(xml_out, "hsum", hsum);
  pop(xml_out);

  multi1d<LatticeReal> meson_phases(Nd);

  MesPhas(meson_phases, 3);
  // Time to bolt
  Chroma::finalize();
  exit(0);

}
